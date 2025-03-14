/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kinematicSingleLayer.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcVolumeIntegrate.H"
#include "fvcFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "mappedWallPolyPatch.H"
#include "mapDistribute.H"
#include "filmThermoModel.H"

// #include "cachedRandom.H" // kvm
#include "normal.H" // kvm
#include "mathematicalConstants.H" // kvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kinematicSingleLayer, 0);

addToRunTimeSelectionTable(surfaceFilmRegionModel, kinematicSingleLayer, mesh);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool kinematicSingleLayer::read()
{
    if (surfaceFilmRegionModel::read())
    {
        const dictionary& solution = this->solution().subDict("PISO");
        solution.readEntry("momentumPredictor", momentumPredictor_);
        solution.readIfPresent("nOuterCorr", nOuterCorr_);
        solution.readEntry("nCorr", nCorr_);
        solution.readEntry("nNonOrthCorr", nNonOrthCorr_);

        return true;
    }

    return false;
}


void kinematicSingleLayer::correctThermoFields()
{
    rho_ == filmThermo_->rho();
    mu_ == filmThermo_->mu();
    sigma_ == filmThermo_->sigma();
}


// kvm
void kinematicSingleLayer::tabAdd(){
    fmtab+="\t";
    return;
}

void kinematicSingleLayer::tabSubtract(){
    fmtab.erase(fmtab.length()-1,fmtab.length());
    return;
}

void kinematicSingleLayer::resetPrimaryRegionSourceTerms()
{
    DebugInFunction << endl;

    rhoSpPrimary_ == dimensionedScalar(rhoSp_.dimensions(), Zero);
    USpPrimary_ == dimensionedVector(USp_.dimensions(), Zero);
    pSpPrimary_ == dimensionedScalar(pSp_.dimensions(), Zero);
}


void kinematicSingleLayer::transferPrimaryRegionThermoFields()
{
    DebugInFunction << endl;

    // Update fields from primary region via direct mapped
    // (coupled) boundary conditions
    UPrimary_.correctBoundaryConditions();
    pPrimary_.correctBoundaryConditions();
    rhoPrimary_.correctBoundaryConditions();
    muPrimary_.correctBoundaryConditions();
}


void kinematicSingleLayer::transferPrimaryRegionSourceFields()
{
    DebugInFunction << endl;

    // Alex
    volScalarField::Boundary& impMassPrimaryBf =
        impMassPrimary_.boundaryFieldRef();
    volVectorField::Boundary& impMomentumPrimaryBf =
        impMomentumPrimary_.boundaryFieldRef();
    volScalarField::Boundary& impPressurePrimaryBf =
        impPressurePrimary_.boundaryFieldRef();
    // Alex

    volScalarField::Boundary& rhoSpPrimaryBf =
        rhoSpPrimary_.boundaryFieldRef();

    volVectorField::Boundary& USpPrimaryBf =
        USpPrimary_.boundaryFieldRef();

    volScalarField::Boundary& pSpPrimaryBf =
        pSpPrimary_.boundaryFieldRef();

    // Convert accumulated source terms into per unit area per unit time
    const scalar deltaT = time_.deltaTValue();
    forAll(rhoSpPrimary_.boundaryField(), patchi)
    {
        scalarField rpriMagSfdeltaT
        (
            (1.0/deltaT)
           /primaryMesh().magSf().boundaryField()[patchi]
        );

        // Alex
        rhoSpPrimaryBf[patchi] = impRelaxFactor_*impMassPrimaryBf[patchi]*rpriMagSfdeltaT;
        impMassPrimaryBf[patchi] *= (1.0-impRelaxFactor_);

        USpPrimaryBf[patchi] = impRelaxFactor_*impMomentumPrimaryBf[patchi]*rpriMagSfdeltaT;
        impMomentumPrimaryBf[patchi] *= (1.0-impRelaxFactor_);

        pSpPrimaryBf[patchi] = impRelaxFactor_*impPressurePrimaryBf[patchi]*rpriMagSfdeltaT;
        impPressurePrimaryBf[patchi] *= (1.0-impRelaxFactor_);
        // Alex
    }

    // Retrieve the source fields from the primary region via direct mapped
    // (coupled) boundary conditions
    // - fields require transfer of values for both patch AND to push the
    //   values into the first layer of internal cells
    rhoSp_.correctBoundaryConditions();
    USp_.correctBoundaryConditions();
    pSp_.correctBoundaryConditions();

    // update addedMassTotal counter
    if (time().writeTime())
    {
        if (debug)
        {
            rhoSp_.write();
            USp_.write();
            pSp_.write();
        }

        scalar addedMassTotal = 0.0;
        outputProperties().readIfPresent("addedMassTotal", addedMassTotal);
        addedMassTotal += returnReduce(addedMassTotal_, sumOp<scalar>());
        outputProperties().add("addedMassTotal", addedMassTotal, true);
        addedMassTotal_ = 0.0;
    }
}

// Alex
// Get information from coupled pyrolysis patch (melting)
void kinematicSingleLayer::transferPyrolysisRegionSourceFields()
{
    DebugInFunction << endl;

    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchi];

        //Info << "Getting melt from patch: " << pp.name() << endl;

        // Copy boundary field to internal cells
        UIndirectList<scalar>(rhoSpPyr_, pp.faceCells()) = rhoSpPyr_.boundaryField()[patchi];

        //Info << rhoSpPyr_ << endl;
    }

    rhoSpPyr_.correctBoundaryConditions();

    return;
}

tmp<volScalarField> kinematicSingleLayer::pu()
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":pu",
                time_.timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pPrimary_                  // pressure (mapped from primary region)
          - pSp_                           // accumulated particle impingement
          - fvc::laplacian(sigma_, delta_) // surface tension
        )
    );
}


tmp<volScalarField> kinematicSingleLayer::pp()
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":pp",
                time_.timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -rho_*gNormClipped() // hydrostatic effect only
        )
    );
}


void kinematicSingleLayer::correctAlpha()
{
    alpha_ == pos(delta_ - deltaSmall_);
}


void kinematicSingleLayer::updateSubmodels()
{
    DebugInFunction << endl;

    // Update injection model - mass returned is mass available for injection
    injection_.correct(availableMass_, cloudMassTrans_, cloudDiameterTrans_);

    // Update transfer model - mass returned is mass available for transfer
    transfer_.correct(availableMass_, cloudMassTrans_);

    // Alex
    // For diagnostics, pull out seperately
    rhoSpDrip_ = cloudMassTrans_/magSf()/time().deltaT();
    rhoSpDrip_.correctBoundaryConditions();

    // Update mass source field
    rhoSp_ += cloudMassTrans_/magSf()/time().deltaT();

    turbulence_->correct();
}


void kinematicSingleLayer::continuityCheck()
{
    const volScalarField deltaRho0(deltaRho_);

    solveContinuity();

    if (debug)
    {
        const volScalarField mass(deltaRho_*magSf());
        const dimensionedScalar totalMass =
            fvc::domainIntegrate(mass)
          + dimensionedScalar("SMALL", dimMass*dimVolume, ROOTVSMALL);

        const scalar sumLocalContErr =
            (
                fvc::domainIntegrate(mag(mass - magSf()*deltaRho0))/totalMass
            ).value();

        const scalar globalContErr =
            (
                fvc::domainIntegrate(mass - magSf()*deltaRho0)/totalMass
            ).value();

        cumulativeContErr_ += globalContErr;

        InfoInFunction
            << "Surface film: " << type() << nl
            << "    time step continuity errors: sum local = "
            << sumLocalContErr << ", global = " << globalContErr
            << ", cumulative = " << cumulativeContErr_ << endl;
    }
}


void kinematicSingleLayer::solveContinuity()
{
    DebugInFunction << endl;

    solve
    (
        fvm::ddt(deltaRho_)
      + fvc::div(phi_)
     ==
      - rhoSp_    // Primary region source term
      + rhoSpPyr_ // Pyrolysis region source term (melt)
    );
}


void kinematicSingleLayer::updateSurfaceVelocities()
{
    // Push boundary film velocity values into internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchi];
        UIndirectList<vector>(Uw_, pp.faceCells()) =
            U_.boundaryField()[patchi];
    }
    Uw_ -= nHat()*(Uw_ & nHat());
    Uw_.correctBoundaryConditions();

    // apply quadratic profile to surface velocity // kvm
    Us_ = 2.0*U_;
    // Limit velocity to reasonable values // kvm
    const scalar Umax=10.0;
    label limitCount=0;
    forAll(Us_,cellI){
        for(label i=0;i<3;i++){
            if(Us_[cellI][i] > Umax){
                Us_[cellI][i] = Umax;
                limitCount++;
            }
            if(Us_[cellI][i] < -Umax){
                Us_[cellI][i] = -Umax;
                limitCount++;
            }
        }
    }
    reduce(limitCount,sumOp<label>());
    if(limitCount>0){
        Info << "limiting " << limitCount << " velocity components" << nl;
    }

    Us_.correctBoundaryConditions();
}


tmp<Foam::fvVectorMatrix> kinematicSingleLayer::solveMomentum
(
    const volScalarField& pu,
    const volScalarField& pp
)
{
    DebugInFunction << endl;
    updateSurfaceVelocities(); // kvm

    // Momentum
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(deltaRho_, U_)
      + fvm::div(phi_, U_)
     ==
      - USp_
   // - fvm::SuSp(rhoSp_, U_)
      - rhoSp_*U_
      + forces_.correct(U_)
      + turbulence_->Su(U_)
    );

    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    if (momentumPredictor_)
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
              - fvc::interpolate(delta_)
              * (
                    regionMesh().magSf()
                  * (
                        fvc::snGrad(pu, "snGrad(p)")
                      + fvc::snGrad(pp, "snGrad(p)")*fvc::interpolate(delta_)
                      + fvc::snGrad(delta_)*fvc::interpolate(pp)
                    )
                  - fvc::flux(rho_*gTan())
                )
            )
        );

        // Remove any patch-normal components of velocity
        U_ -= nHat()*(nHat() & U_);
        U_.correctBoundaryConditions();
    }

    return tUEqn;
}


void kinematicSingleLayer::solveThickness
(
    const volScalarField& pu,
    const volScalarField& pp,
    fvVectorMatrix& UEqn
)
{
    DebugInFunction << endl;

    volScalarField rUA(1.0/UEqn.A());
    U_ = rUA*UEqn.H();

    surfaceScalarField deltarUAf(fvc::interpolate(delta_*rUA));
    surfaceScalarField rhof(fvc::interpolate(rho_));

    surfaceScalarField phiAdd
    (
        "phiAdd",
        regionMesh().magSf()
      * (
            fvc::snGrad(pu, "snGrad(p)")
          + fvc::snGrad(pp, "snGrad(p)")*fvc::interpolate(delta_)
        )
      - fvc::flux(rho_*gTan())
    );
    constrainFilmField(phiAdd, 0.0);

    surfaceScalarField phid
    (
        "phid",
        fvc::flux(U_*rho_) - deltarUAf*phiAdd*rhof
    );
    constrainFilmField(phid, 0.0);

    surfaceScalarField ddrhorUAppf
    (
        "deltaCoeff",
        fvc::interpolate(delta_)*deltarUAf*rhof*fvc::interpolate(pp)
    );

    regionMesh().setFluxRequired(delta_.name());

    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        // Film thickness equation
        fvScalarMatrix deltaEqn
        (
            fvm::ddt(rho_, delta_)
          + fvm::div(phid, delta_)
          - fvm::laplacian(ddrhorUAppf, delta_)
         ==
          - rhoSp_   // Primary region source terms
          + rhoSpPyr_ // Pyrolysis region source term (Melt)
        );

        deltaEqn.solve();

        if (nonOrth == nNonOrthCorr_)
        {
            phiAdd +=
                fvc::interpolate(pp)
              * fvc::snGrad(delta_)
              * regionMesh().magSf();

            phi_ == deltaEqn.flux();
        }
    }

    // Bound film thickness by a minimum of zero
    delta_.clamp_min(0);

    // Update U field
    U_ -= fvc::reconstruct(deltarUAf*phiAdd);

    // Remove any patch-normal components of velocity
    U_ -= nHat()*(nHat() & U_);

    U_.correctBoundaryConditions();

    // Continuity check
    continuityCheck();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kinematicSingleLayer::kinematicSingleLayer
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    surfaceFilmRegionModel(modelType, mesh, g, regionType),

    momentumPredictor_(solution().subDict("PISO").lookup("momentumPredictor")),
    nOuterCorr_(solution().subDict("PISO").getOrDefault("nOuterCorr", 1)),
    nCorr_(solution().subDict("PISO").get<label>("nCorr")),
    nNonOrthCorr_
    (
        solution().subDict("PISO").get<label>("nNonOrthCorr")
    ),

    cumulativeContErr_(0.0),

    deltaSmall_("deltaSmall", dimLength, SMALL),
    deltaCoLimit_(solution().getOrDefault<scalar>("deltaCoLimit", 1e-4)),
    impRelaxFactor_(coeffs_.lookupOrDefault("impingementRelaxationFactor",1.0)), // Alex
    // Ning
    seamArrayX_(coeffs_.lookupOrDefault("seamArrayX",(0.0))),
    seamArrayY_(coeffs_.lookupOrDefault("seamArrayY",(0.0))),
    rackTiers_(coeffs_.lookupOrDefault("rackTiers", 3)),
    rackElevated_(coeffs_.lookupOrDefault("rackElevated", 8)),
    rackTierSeperationDistance_(coeffs_.lookupOrDefault("rackTierSeperationDistance", 60)),

    seamCells_
    (
        IOobject
        (
            "seamCells",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    // Ning
    rho_
    (
        IOobject
        (
            "rhof",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimDensity, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    mu_
    (
        IOobject
        (
            "muf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure*dimTime, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    sigma_
    (
        IOobject
        (
            "sigmaf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass/sqr(dimTime), Zero),
        fvPatchFieldBase::zeroGradientType()
    ),

    delta_
    (
        IOobject
        (
            "deltaf",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    alpha_
    (
        IOobject
        (
            "alpha",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    U_
    (
        IOobject
        (
            "Uf",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    Us_
    (
        IOobject
        (
            "Usf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_,
        fvPatchFieldBase::zeroGradientType()
    ),
    Uw_
    (
        IOobject
        (
            "Uwf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_,
        fvPatchFieldBase::zeroGradientType()
    ),
    deltaRho_
    (
        IOobject
        (
            delta_.name() + "*" + rho_.name(),
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(delta_.dimensions()*rho_.dimensions(), Zero),
        fvPatchFieldBase::zeroGradientType()
    ),

    phi_
    (
        IOobject
        (
            "phi",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT, // kvm
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimLength*dimMass/dimTime, Zero)
    ),

    primaryMassTrans_
    (
        IOobject
        (
            "primaryMassTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    cloudMassTrans_
    (
        IOobject
        (
            "cloudMassTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    cloudDiameterTrans_
    (
        IOobject
        (
            "cloudDiameterTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("minus1", dimLength, -1.0),
        fvPatchFieldBase::zeroGradientType()
    ),
    USp_
    (
        IOobject
        (
            "USpf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedVector(dimMass*dimVelocity/dimArea/dimTime, Zero),
        this->mappedPushedFieldPatchTypes<vector>()
    ),
    pSp_
    (
        IOobject
        (
            "pSpf",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, Zero),
        this->mappedPushedFieldPatchTypes<scalar>()
    ),
    rhoSp_
    (
        IOobject
        (
            "rhoSpf",
            time_.timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass/dimTime/dimArea, Zero),
        this->mappedPushedFieldPatchTypes<scalar>()
    ),

    USpPrimary_
    (
        IOobject
        (
            USp_.name(), // must have same name as USp_ to enable mapping
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedVector(USp_.dimensions(), Zero)
    ),
    pSpPrimary_
    (
        IOobject
        (
            pSp_.name(), // must have same name as pSp_ to enable mapping
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(pSp_.dimensions(), Zero)
    ),
    rhoSpPrimary_
    (
        IOobject
        (
            rhoSp_.name(), // must have same name as rhoSp_ to enable mapping
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(rhoSp_.dimensions(), Zero)
    ),
    // Alex
    impMassPrimary_
    (
        IOobject
        (
            "impMass", 
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(dimMass, Zero)
    ),
    impMomentumPrimary_
    (
        IOobject
        (
            "impMomentum", 
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedVector(dimMass*dimVelocity, Zero)
    ),
    impPressurePrimary_
    (
        IOobject
        (
            "impPressure", 
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    // Mass flux due to gasification (phaseChange submodel only)
    rhoSpGas_
    (
        IOobject
        (
            "filmGasMassFlux", 
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimMass/dimArea/dimTime, Zero)
    ),
    // Mass flux from pyrolysis region (melting)
    rhoSpPyr_
    (
        IOobject
        (
            "rhoSpPyr",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(rhoSp_.dimensions(), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    // Mass flux droplet impingement
    rhoSpImp_
    (
        IOobject
        (
            "rhoSpImp",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(rhoSp_.dimensions(), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    // Mass flux droplet dripping
    rhoSpDrip_
    (
        IOobject
        (
            "rhoSpDrip",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(rhoSp_.dimensions(), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    // Alex
    UPrimary_
    (
        IOobject
        (
            "U", // must have same name as U to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedVector(dimVelocity, Zero),
        this->mappedFieldAndInternalPatchTypes<vector>()
    ),
    pPrimary_
    (
        IOobject
        (
            "p", // must have same name as p to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, Zero),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),
    rhoPrimary_
    (
        IOobject
        (
            "rho", // must have same name as rho to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimDensity, Zero),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),
    muPrimary_
    (
        IOobject
        (
            "thermo:mu", // must have same name as mu to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure*dimTime, Zero),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),

    filmThermo_(filmThermoModel::New(*this, coeffs_)),

    availableMass_(regionMesh().nCells(), Zero),

    injection_(*this, coeffs_),

    transfer_(*this, coeffs_),

    turbulence_(filmTurbulenceModel::New(*this, coeffs_)),

    forces_(*this, coeffs_),

    pyrolysisMelt_(coeffs().lookupOrDefault("pyrolysisMelt",false)), // Alex

    addedMassTotal_(0.0)
{
    if (readFields)
    {
        transferPrimaryRegionThermoFields();

        correctAlpha();

        correctThermoFields();

        // delta needs to be initialized before deltaRho_ calculation (call bc) // kvm
        deltaRho_ == delta_*rho_;
        // U_ needs to be initialized before phi_ calculation (call bc) // kvm
        surfaceScalarField phi0
        (
            IOobject
            (
                "phi",
                time().timeName(),
                regionMesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            fvc::flux(deltaRho_*U_)
        );

        phi_ == phi0;
    }
    forAll(delta_,cellI)
    {
        scalar cellX(regionMesh().C()[cellI].x());
        scalar cellY(regionMesh().C()[cellI].y());
        scalar cellZ(regionMesh().C()[cellI].z());
        for(int i=0;i<rackTiers_;i++)
        {
            scalar seamZ(0.0254*(rackElevated_ + 47.0 + rackTierSeperationDistance_*i));
            if(fabs(cellZ - seamZ) < 0.01)
            {
                forAll(seamArrayX_, iX)
                {
                    if(fabs(cellX - seamArrayX_[iX]) < 0.025) seamCells_[cellI] = 1.0;
                }
                forAll(seamArrayY_, iY)
                {
                    if(fabs(cellY - seamArrayY_[iY]) < 0.025) seamCells_[cellI] = 1.0;
                }
            }
        }
    }
    seamCells_.correctBoundaryConditions();


    impRelaxFactor_ = min(1.0,max(0.0,impRelaxFactor_));
    Info << "\tFilm model relaxation factor for impingement sources: " << impRelaxFactor_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

kinematicSingleLayer::~kinematicSingleLayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void kinematicSingleLayer::addSources
(
    const label patchi,
    const label facei,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    DebugInFunction
        << "\nSurface film: " << type() << ": adding to film source:" << nl
        << "    mass     = " << massSource << nl
        << "    momentum = " << momentumSource << nl
        << "    pressure = " << pressureSource << endl;

    // Alex
    impMassPrimary_.boundaryFieldRef()[patchi][facei] -= massSource;
    impMomentumPrimary_.boundaryFieldRef()[patchi][facei] -= momentumSource;
    impPressurePrimary_.boundaryFieldRef()[patchi][facei] -= pressureSource;

    addedMassTotal_ += massSource;
}


void kinematicSingleLayer::preEvolveRegion()
{
    DebugInFunction << endl;

    surfaceFilmRegionModel::preEvolveRegion();

    transferPrimaryRegionThermoFields();

    correctThermoFields();

    transferPrimaryRegionSourceFields();

    // Alex
    if (pyrolysisMelt_)
    {
        Info << "Transfer melt from coupled pyrolysis region" << endl;
        transferPyrolysisRegionSourceFields();
    }

    //- At this stage rhoSp_ has only picked up impingement information 
    //- from transferPrimaryRegionSourceFields()
    //- Later, dripping and other stuff will be lumped in via updateSubmodels()
    //- So, the idea here is to cache the impingement information
    rhoSpImp_ = rhoSp_;
    rhoSpImp_.correctBoundaryConditions();

    // Reset transfer fields
    //availableMass_ = mass();
    availableMass_ = netMass();
    cloudMassTrans_ == dimensionedScalar(dimMass, Zero);
    cloudDiameterTrans_ == dimensionedScalar(dimLength, -1.0); // kvm
    primaryMassTrans_ == dimensionedScalar(dimMass, Zero);
    // Alex
}


void kinematicSingleLayer::evolveRegion()
{
    DebugInFunction << endl;

    // Update film coverage indicator
    correctAlpha();

    // Update sub-models to provide updated source contributions
    updateSubmodels();

    // Solve continuity for deltaRho_
    solveContinuity();

    // Implicit pressure source coefficient - constant
    tmp<volScalarField> tpp(this->pp());

    for (int oCorr=1; oCorr<=nOuterCorr_; oCorr++)
    {
        // Explicit pressure source contribution - varies with delta_
        tmp<volScalarField> tpu(this->pu());

        // Solve for momentum for U_
        tmp<fvVectorMatrix> tUEqn = solveMomentum(tpu(), tpp());
        fvVectorMatrix& UEqn = tUEqn.ref();

        // Film thickness correction loop
        for (int corr=1; corr<=nCorr_; corr++)
        {
            // Solve thickness for delta_
            solveThickness(tpu(), tpp(), UEqn);
        }
    }

    // Update deltaRho_ with new delta_
    deltaRho_ == delta_*rho_;

    // Update film wall and surface velocities // kvm
    updateSurfaceVelocities(); // kvm

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();
}

void kinematicSingleLayer::postEvolveRegion()
{
    DebugInFunction << endl;

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();
}


scalar kinematicSingleLayer::CourantNumber() const
{
    scalar CoNum = 0.0;

    if (regionMesh().nInternalFaces() > 0)
    {
        const scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi_))().primitiveField()
          / (deltaRho_.primitiveField() + ROOTVSMALL)
        );

        forAll(delta_, i)
        {
            if ((delta_[i] > deltaCoLimit_) && (alpha_[i] > 0.5))
            {
                CoNum = max(CoNum, sumPhi[i]/(delta_[i]*magSf()[i]));
            }
        }

        CoNum *= 0.5*time_.deltaTValue();
    }

    reduce(CoNum, maxOp<scalar>());

    Info<< "Film max Courant number: " << CoNum << endl;

    return CoNum;
}


const volVectorField& kinematicSingleLayer::U() const
{
    return U_;
}


const volVectorField& kinematicSingleLayer::Us() const
{
    return Us_;
}


const volVectorField& kinematicSingleLayer::Uw() const
{
    return Uw_;
}


const volScalarField& kinematicSingleLayer::deltaRho() const
{
    return deltaRho_;
}


const surfaceScalarField& kinematicSingleLayer::phi() const
{
    return phi_;
}


const volScalarField& kinematicSingleLayer::rho() const
{
    return rho_;
}


const volScalarField& kinematicSingleLayer::T() const
{
    FatalErrorInFunction
        << "T field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& kinematicSingleLayer::Ts() const
{
    FatalErrorInFunction
        << "Ts field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& kinematicSingleLayer::Tw() const
{
    FatalErrorInFunction
        << "Tw field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& kinematicSingleLayer::hs() const
{
    FatalErrorInFunction
        << "hs field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& kinematicSingleLayer::Cp() const
{
    FatalErrorInFunction
        << "Cp field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const volScalarField& kinematicSingleLayer::kappa() const
{
    FatalErrorInFunction
        << "kappa field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


tmp<volScalarField> kinematicSingleLayer::primaryMassTrans() const
{
    return primaryMassTrans_;
}


const volScalarField& kinematicSingleLayer::cloudMassTrans() const
{
    return cloudMassTrans_;
}


const volScalarField& kinematicSingleLayer::cloudDiameterTrans() const
{
    return cloudDiameterTrans_;
}


void kinematicSingleLayer::info()
{
    Info<< "\nSurface film: " << type() << endl;

    const scalarField& deltaInternal = delta_.primitiveField(); // kvm
    const vectorField& Uinternal = U_.primitiveField(); // kvm
    scalar addedMassTotal = 0.0;
    outputProperties().readIfPresent("addedMassTotal", addedMassTotal);
    addedMassTotal += returnReduce(addedMassTotal_, sumOp<scalar>());

    Info<< indent << "added mass         = " << addedMassTotal << nl
        << indent << "current mass       = "
        << gSum((deltaRho_*magSf())()) << nl
        << indent << "min/max(mag(U))    = " << gMin(mag(Uinternal)) << ", "
        << gMax(mag(Uinternal)) << nl
        << indent << "min/max(delta)     = " << gMin(deltaInternal) << ", "
        << gMax(deltaInternal) << nl
        << indent << "coverage           = "
        << gSum(alpha_.primitiveField()*magSf())/gSum(magSf()) <<  nl;

    injection_.info(Info);
    transfer_.info(Info);
}


tmp<volScalarField::Internal> kinematicSingleLayer::Srho() const
{
    return tmp<volScalarField::Internal>::New
    (
        IOobject
        (
            typeName + ":Srho",
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        primaryMesh(),
        dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
    );
}


tmp<volScalarField::Internal> kinematicSingleLayer::Srho
(
    const label i
) const
{
    return tmp<volScalarField::Internal>::New
    (
        IOobject
        (
            typeName + ":Srho(" + Foam::name(i) + ")",
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        primaryMesh(),
        dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
    );
}


tmp<volScalarField::Internal> kinematicSingleLayer::Sh() const
{
    return tmp<volScalarField::Internal>::New
    (
        IOobject
        (
            typeName + ":Sh",
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        primaryMesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
