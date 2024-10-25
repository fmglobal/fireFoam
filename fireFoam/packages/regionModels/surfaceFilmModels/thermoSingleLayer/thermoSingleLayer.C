/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "thermoSingleLayer.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcFlux.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "mappedFieldFvPatchField.H"
#include "mapDistribute.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// Sub-models
#include "filmThermoModel.H"
#include "filmViscosityModel.H"
#include "heatTransferModel.H"
#include "phaseChangeModel.H"
#include "massAbsorptionModel.H"
#include "filmRadiationModel.H"
#include "pyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermoSingleLayer, 0);

addToRunTimeSelectionTable(surfaceFilmRegionModel, thermoSingleLayer, mesh);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

wordList thermoSingleLayer::hsBoundaryTypes()
{
    wordList bTypes(T_.boundaryField().types());
    forAll(bTypes, patchi)
    {
        if
        (
            T_.boundaryField()[patchi].fixesValue()
         || isA<mixedFvPatchScalarField>(T_.boundaryField()[patchi])
         || isA<mappedFieldFvPatchField<scalar>>(T_.boundaryField()[patchi])
        )
        {
            bTypes[patchi] = fixedValueFvPatchField<scalar>::typeName;
        }
    }

    return bTypes;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool thermoSingleLayer::read()
{
    // No additional properties to read
    return kinematicSingleLayer::read();
}


void thermoSingleLayer::resetPrimaryRegionSourceTerms()
{
    DebugInFunction << endl;

    kinematicSingleLayer::resetPrimaryRegionSourceTerms();

    hsSpPrimary_ == dimensionedScalar(hsSp_.dimensions(), Zero);
}


void thermoSingleLayer::correctThermoFields()
{
    // limit temperature to reasonable bounds
    T_.clamp_range(Tbounds_);
    Ts_.clamp_range(Tbounds_);

    rho_ == filmThermo_->rho();
    mu_ == filmThermo_->mu(); // OpenFOAM version missing this, bug, kvm
    sigma_ == filmThermo_->sigma();
    Cp_ == filmThermo_->Cp();
    kappa_ == filmThermo_->kappa();

    rho_.correctBoundaryConditions(); // kvm
    mu_.correctBoundaryConditions(); // kvm
    sigma_.correctBoundaryConditions(); // kvm
    Cp_.correctBoundaryConditions(); // kvm
    kappa_.correctBoundaryConditions(); // kvm

}


void thermoSingleLayer::correctHsForMappedT()
{
    T_.correctBoundaryConditions();

    volScalarField::Boundary& hsBf = hs_.boundaryFieldRef();

    forAll(hsBf, patchi)
    {
        const fvPatchField<scalar>& Tp = T_.boundaryField()[patchi];
        const fvPatchField<scalar>& Pp = this->pPrimary().boundaryField()[patchi]; // Alex 
        if (isA<mappedFieldFvPatchField<scalar>>(Tp))
        {
            // Alex
            forAll(Tp,faceI)
            {
                hsBf[patchi][faceI] = filmThermo_->hs(Pp[faceI],Tp[faceI],true); 
            }
            // Alex
        }
    }
}


void thermoSingleLayer::updateSurfaceTemperatures()
{
    correctHsForMappedT();

    // FM - pull char fraction from pyrolysis region
    charFrac_.correctBoundaryConditions();

    // Push boundary film temperature into wall temperature internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = regionMesh().boundaryMesh()[patchi];
        if(pyrCoupled_){
        // deprecated    // get pyrolysis internal temperature 
        // deprecated    // and put it on coupled patch boundary of film T_ 

        // deprecated     // UIndirectList<scalar>(Tw_, pp.faceCells()) =
        // deprecated     // 		//T_.boundaryFieldRef()[patchi];
        // deprecated     // 		pyrTemperaturePtr_->boundaryFieldRef()[patchi];
        // deprecated     // TODO: how do I do this with the new coupling?

        // deprecated     scalarList Tpyr(pp.faceCells().size(), 0.0);


        // deprecated     typedef regionModels::pyrolysisModels::pyrolysisModel
        // deprecated         pyrolysisModelType;

        // deprecated     const regionModels::regionModel& pyrolysisRegion =
        // deprecated     db().time().lookupObject<regionModels::regionModel>
        // deprecated         (
        // deprecated             "pyrolysisProperties"
        // deprecated             );

        // deprecated     const pyrolysisModelType& pyrolysisModel =
        // deprecated         dynamic_cast<const pyrolysisModelType&>(pyrolysisRegion);
        // deprecated     
        // deprecated     pyrolysisModelType& pyrolysis =
        // deprecated         const_cast<pyrolysisModelType&>(pyrolysisModel);
        // deprecated     
        // deprecated     // internal cell tempertaure must be used for stability
        // deprecated     Tpyr = 
        // deprecated         mapRegionPatchInternalField<scalar>
        // deprecated         (
        // deprecated             pyrolysis,
        // deprecated             "T",
        // deprecated             patchi,
        // deprecated             true
        // deprecated             );
        // deprecated     
        // deprecated     UIndirectList<scalar>(Tw_,pp.faceCells()) = 
        // deprecated         Tpyr;

        // deprecated     // get char fraction from pyrolysis model
        // deprecated     scalarList charFrac(pp.faceCells().size(), 0.0);
        // deprecated     charFrac = 
        // deprecated         mapRegionPatchInternalField<scalar>
        // deprecated         (
        // deprecated             pyrolysis,
        // deprecated             "charFrac",
        // deprecated             patchi,
        // deprecated             true
        // deprecated             );

        // deprecated     UIndirectList<scalar>(charFrac_,pp.faceCells()) = 
        // deprecated         charFrac;

        // deprecated     
        // deprecated }
        // else if(1){
        //     // compute Tw based on 0D lumped capacitance model
        //     forAll(Tw_,i){
        //         Tw_[i]=qFilmToWall_[i]*time_.deltaTValue()/2702.0/949.0/0.0012 + Tw_[i];
        //     }
        //     Info << "max Tw " << tab << db().time().timeName() << tab << gMax(Tw_) << endl;
        }
        else{
            UIndirectList<scalar>(Tw_, pp.faceCells()) =
            		T_.boundaryField()[patchi];
            // set char fraction internal field equal to boundary values
            UIndirectList<scalar>(charFrac_, pp.faceCells()) =
            		charFrac_.boundaryField()[patchi];
        }
    }
    Tw_.correctBoundaryConditions();

    // Update film surface temperature
    Ts_ = T_;
    Ts_.correctBoundaryConditions();
}


tmp<volScalarField> thermoSingleLayer::invertEnthalpy() const
{
    /*
     Newton method used to invert T from hs

     f: true enthalpy - enthalpy calculated from T
     f`: specific heat 
     Uses current film.T_ value for initial guess and film.hs_ for
     true value of enthalpy
     tTn passed out satisfies optimization criteria 
    */

    tmp<volScalarField> tTn
    (
        new volScalarField
        (
            IOobject
            (
                type() + ":Tn",
                time().timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            regionMesh(),
            dimensionedScalar(dimTemperature,Zero),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    scalarField   To = T_.primitiveField();
    scalarField&  Tn = tTn.ref().primitiveFieldRef();
    const scalarField&  h  = hs_.primitiveField();
    const volScalarField& p = this->pPrimary();

    // Conversion criteria
    const scalar rtol(1e-4);
    const label maxIter(20);

    // Newton loop
    forAll(Tn,celli)
    {
        label newtonCount = 0;
        scalar f(Zero);
        scalar residual(Zero);
        do
        {
        // difference between hs solution and hs implied from T_ 

            //Info << "To: " << To[celli] << endl;
            //Info << "hs: " << h[celli] << endl;
            //Info << "hs(p,To): " << filmThermo_->hs(p[celli],Tn[celli]) << endl;
            //Info << "Cp(To): " << filmThermo_->Cp(p[celli],To[celli]) << endl;
            //Info << "hinv: " << hinv << endl;

            // Sensible enthalpy with respect to reference temperature
            const scalar hinv = filmThermo_->hs(p[celli],To[celli],true);
           
            f = - h[celli] + hinv;

            Tn[celli] = To[celli] - f/(filmThermo_->Cp(p[celli],To[celli]));
            
            residual = mag(Tn[celli]-To[celli])/Tn[celli]; 

            ++newtonCount;

            To[celli] = Tn[celli];            

            //Info << "f: " << f << endl;
            //Info << "Tn: " << Tn[celli] << endl;
            //Info << "residual: " << residual << endl;
            //Info << endl;

        }while(mag(residual)>rtol && newtonCount< maxIter);

    }

    Tn.clamp_range(Tbounds_);
    return tTn;
}


void thermoSingleLayer::transferPrimaryRegionThermoFields()
{
    DebugInFunction << endl;

    kinematicSingleLayer::transferPrimaryRegionThermoFields();

    // Update primary region fields on local region via direct mapped (coupled)
    // boundary conditions
    TPrimary_.correctBoundaryConditions();
    forAll(YPrimary_, i)
    {
        YPrimary_[i].correctBoundaryConditions();
    }
}


void thermoSingleLayer::transferPrimaryRegionSourceFields()
{
    DebugInFunction << endl;

    kinematicSingleLayer::transferPrimaryRegionSourceFields();

    volScalarField::Boundary& impEnergyPrimaryBf =
        impEnergyPrimary_.boundaryFieldRef();

    volScalarField::Boundary& hsSpPrimaryBf =
        hsSpPrimary_.boundaryFieldRef();

    // Convert accumulated source terms into per unit area per unit time
    const scalar deltaT = time_.deltaTValue();
    forAll(hsSpPrimaryBf, patchi)
    {
        scalarField rpriMagSfdeltaT
        (
            (1.0/deltaT)/primaryMesh().magSf().boundaryField()[patchi]
        );

        hsSpPrimaryBf[patchi] = impRelaxFactor_*impEnergyPrimaryBf[patchi]*rpriMagSfdeltaT; // Alex
        impEnergyPrimaryBf[patchi] *= (1-impRelaxFactor_); // Alex
    }

    // Retrieve the source fields from the primary region via direct mapped
    // (coupled) boundary conditions
    // - fields require transfer of values for both patch AND to push the
    //   values into the first layer of internal cells
    hsSp_.correctBoundaryConditions();

    // Apply enthalpy source as difference between incoming and actual states
    hsSp_ -= rhoSp_*hs_; 
}


// FM - needed for melt transfer
void thermoSingleLayer::transferPyrolysisRegionSourceFields()
{
    // rhoSpPyr_ updated here
    kinematicSingleLayer::transferPyrolysisRegionSourceFields();

    // The melt comes in at the melt temperature, which is same as the wall temperature
    // mapped from the pyrolysis region
    tmp<volScalarField> hmelt(filmThermo_->hs(Tw_,true)); // - filmThermo_->hs(T_,false));
    hsSpPyr_ = rhoSpPyr_*hmelt;
    hsSpPyr_.correctBoundaryConditions();

    return;
}


void thermoSingleLayer::correctAlpha()
{
    if (hydrophilic_)
    {
        const scalar hydrophilicDry = hydrophilicDryScale_*deltaWet_;
        const scalar hydrophilicWet = hydrophilicWetScale_*deltaWet_;

        forAll(alpha_, i)
        {
            if ((alpha_[i] < 0.5) && (delta_[i] > hydrophilicWet))
            {
                alpha_[i] = 1.0;
            }
            else if ((alpha_[i] > 0.5) && (delta_[i] < hydrophilicDry))
            {
                alpha_[i] = 0.0;
            }
        }

        alpha_.correctBoundaryConditions();
    }
    else
    {
        alpha_ ==
            pos0(delta_ - dimensionedScalar("deltaWet", dimLength, deltaWet_));
    }
}

// FM 
void thermoSingleLayer::updatePhaseChange() 
{

    // Update vaporization
    phaseChange_->correct
    (
        time_.deltaTValue(),
        availableMass_,
        primaryMassTrans_,
        primaryEnergyTrans_
    );


    rhoSpGas_ = primaryMassTrans_/this->magSf()/time_.deltaT();

}

// FM
void thermoSingleLayer::updateConvectiveFlux()
{
    // Update heat transfer to wall (used in film/pyrolysis coupling)
    // heat flow out of film is positive
    qFilmToWall_ = alpha_*htcw_->h()*(T_ - Tw_);
    qFilmToWall_.correctBoundaryConditions();

    // Update heat transfer from gas phase 
    // heat flow out of film is positive
    qGasToFilm_ = alpha_*htcs_->h()*(T_ - TPrimary_);
    qGasToFilm_.correctBoundaryConditions();

}

void thermoSingleLayer::updateHeatFluxCorrectionZone()
{
/*
    Tags cells for the optional additional of heat flux source terms in
    temperature BC.

    The idea is to find the "pyrolysis front", then tag cells on the unburnt side
    of the front, within some distance, which is then used as a mask for applying
    the additional heat flux term.

    The method is as follows:

    Consider a 1D mesh of 4 cells:

       cell1    cell2    cell3    cell4
    |    x   |    o   |    o   |    o   | 

    The 'x' indicates MLR>MLR_crit in cell1, which is the pyrolysis front

    First we tag cell1 as being in the burnt zone:
       cell1    cell2    cell3    cell4
    |  bz=1  |  bz=0  |  bz=0  |  bz=0  | 

    Then we tag cell2 as being in the correction zone:
       cell1    cell2    cell3    cell4
    |  bz=1  |  bz=0  |  bz=0  |  bz=0  |
    |  cz=0  |  cz=1  |  cz=0  |  cz=0  |

    Then we search for cells next to cells tagged with cz=1
    and not in the burnt zone and tag them with cz=2

    cell1    cell2    cell3    cell4
    |  bz=1  |  bz=0  |  bz=0  |  bz=0  |
    |  cz=0  |  cz=1  |  cz=2  |  cz=0  |

    We repeat this process, tagging values of correctionZone up to
    and including the user-specified value of preHeatingCells_ (here, it's 2)

    Note that the tag value is incremented to track the width of the correction zone
    but that all tagged cells will have same heat flux applied

    The correctionZone_ field is then referenced in the coupled temperature
    BC and where correctionZone_ >= 1, the additional heat flux is applied. 
*/

    phiSolid_.correctBoundaryConditions();

    massLossRate_ = phiSolid_;

    // push values from patch to film internal field
    for (label i=0; i<intCoupledPatchIDs_.size(); i++)
    {
        const label patchi = intCoupledPatchIDs_[i];
        scalarField& phiSolidW = phiSolid_.boundaryFieldRef()[patchi];
        const labelUList& faceCells = regionMesh().boundary()[patchi].faceCells();
        forAll(phiSolidW,facei)
        {
            massLossRate_[faceCells[facei]] = phiSolidW[facei];
        }
    }

    massLossRate_.correctBoundaryConditions();

    const vectorField& Cv = regionMesh().C();
    const scalarField& cellVolume = regionMesh().V();

    forAll(massLossRate_, cellI)
    {
        correctionZone_[cellI] = 0;

        // Tag cells to identify the 'pyrolysis front'
        if(massLossRate_[cellI] > 3.0*burntMlrCrit_)
        {
            if(burntZone_[cellI] == 0)
            {
                burntZone_[cellI] = 1.0;
            }
        }


        // Looking for cells on the un-reacted side of the pyrolysis front
        if(massLossRate_[cellI] < burntMlrCrit_)
        {
            // Neighbour cells to current cell
            const labelList& nbrCL = regionMesh().cellCells()[cellI];
            forAll(nbrCL,nI)
            {
                label nbrCellID = nbrCL[nI];
                // If neighbour cell is reacting and current cell is not in burnt zone
                // , then tag it for correction
                if((massLossRate_[nbrCellID] > burntMlrCrit_) && (burntZone_[cellI] == 0))
                {
                    correctionZone_[cellI] = 1;
                }
            }
        }
    }

    // Check across processor boundaries to apply correction
    forAll(massLossRate_.boundaryField(),pI)
    {
        const fvPatchField<scalar>& mlrB = massLossRate_.boundaryField()[pI];
    
        if(mlrB.type() == "processor")
        {
            const labelUList& faceCells = mlrB.patch().faceCells();
            forAll(mlrB,fi)
            {
                label cellO = faceCells[fi];
                if((burntZone_[cellO] == 0) && (massLossRate_[cellO] < burntMlrCrit_) && (mlrB[fi] > burntMlrCrit_))
                {
                    correctionZone_[cellO] = 1;
                }
            }
        }
    }

    correctionZone_.correctBoundaryConditions();
    burntZone_.correctBoundaryConditions();

    // Thicken the correction zone by preHeatingCells_ number of cells
    for(int ipc = 1; ipc < preHeatingCells_; ipc++)
    {
        Info<<"CorrectionZone No.: "<<ipc+1<<endl;
        forAll(massLossRate_, cellI)
        {
            // If not in the burnt zone or the correction zone
            if((burntZone_[cellI] == 0) && (correctionZone_[cellI] == 0) && (massLossRate_[cellI] < burntMlrCrit_))
            {
                const labelList& nbrCL = regionMesh().cellCells()[cellI];
                // Loop over neighbor cells
                forAll(nbrCL,nI)
                {
                    label nbrCellID = nbrCL[nI];
                    // If a neighbour is in the correction zone (tagged with number of cells away from the burnt zone)
                    // and its correctionZone tag matches the current level, then add this cell to the correction
                    // zone and increment its value
                    if(correctionZone_[nbrCellID] == ipc)
                    {
                        correctionZone_[cellI] = 1 + ipc;
                    }
                }
            }
        }
    }
    correctionZone_.correctBoundaryConditions();
}

void thermoSingleLayer::updateSubmodels()
{
    DebugInFunction << endl;

    // Update heat transfer coefficient sub-models
    htcs_->correct();
    htcw_->correct();

    // FM 
    // Wrapper to set qFilmToWall and qGasToFilm
    updateConvectiveFlux();

    // Wrapper to phaseChange_->correct
    updatePhaseChange();

    //for diagnostics
    primaryMassTrans_.correctBoundaryConditions();
    primaryEnergyTrans_.correctBoundaryConditions();
    charFrac_.correctBoundaryConditions();

    // Update massAbsorption
    massAbsorption_->correct
    (
        time_.deltaTValue(),
        availableMass_,
        massAbs_,
        energyAbs_
    );

    //for diagnostics
    massAbs_.correctBoundaryConditions();
    energyAbs_.correctBoundaryConditions();

    // Update radiation
    radiation_->correct();
    // This is picked up by regionCoupling T BC 
    radTransmitted_ = radiation_->transmitted();
    radTransmitted_.correctBoundaryConditions();

    // This is needed specifically by CUP pyrolysis model 
    // (re-evaluate this decision)
	transmissivity_ = radiation_->transmissivity();
	transmissivity_.correctBoundaryConditions();

    // Kinematic stuff

    // Update injection model - mass returned is mass available for injection
    injection_.correct(availableMass_, cloudMassTrans_, cloudDiameterTrans_);

    // Update transfer model - mass returned is mass available for transfer
    transfer_.correct(availableMass_, cloudMassTrans_);

    // For diagnostics, pull out seperately
    rhoSpDrip_ = cloudMassTrans_/magSf()/time().deltaT();
    rhoSpDrip_.correctBoundaryConditions();

    // Update mass source field
    rhoSp_ += cloudMassTrans_/magSf()/time().deltaT();

    turbulence_->correct();

    // - cloudEnergyTrans 
    hsSpDrip_ = rhoSpDrip_*hs_;

    // Update transfer model - mass returned is mass available for transfer
    transfer_.correct(availableMass_, primaryMassTrans_, primaryEnergyTrans_);

    // Cache phase change balues
    rhoSpGas_ = primaryMassTrans_/magSf()/time().deltaT();
    hsSpGas_ = primaryEnergyTrans_/magSf()/time().deltaT();

    // Update source fields (phase change)
    hsSp_ += primaryEnergyTrans_/magSf()/time().deltaT();
    rhoSp_ += primaryMassTrans_/magSf()/time().deltaT();

    // Update source fields (mass absorption)
    rhoSp_ += massAbs_/magSf()/time().deltaT();

    // Vapour recoil pressure (can become unstable for wild oscillations in vaporization rate)
    // pSp_ -= sqr(primaryMassTrans_/magSf()/time_.deltaT())/2.0/rhoPrimary_;
    // Info << "vaporRecoilPressure " << gMin(pSp_) << " " << gAverage(pSp_) << " " << gMax(pSp_) << nl;
    // FM

    // Enhanced lateral flame spread
    updateHeatFluxCorrectionZone();  
}


tmp<fvScalarMatrix> thermoSingleLayer::q(volScalarField& hs) const
{

    Info << "htcs_->h(): " << tab << db().time().timeName() << tab << max(htcs_->h()) << endl;
    Info << "htcw_->h(): " << tab << db().time().timeName() << tab << max(htcw_->h()) << endl;

    scalar maxTw = max(Tw_.primitiveField());
    reduce(maxTw,maxOp<scalar>());
    Info << "max(Tw_): " << tab << db().time().timeName() << tab << maxTw << endl;

    scalar maxTPrimary = max(TPrimary_.primitiveField());
    reduce(maxTPrimary,maxOp<scalar>());
    Info << "max(TPrimary_): " << tab << db().time().timeName() << tab << maxTPrimary << endl;

    return
    (
        // Heat-transfer to the primary region
      - fvm::Sp(htcs_->h()/Cp_, hs)
      + htcs_->h()*(hs/Cp_ + alpha_*(TPrimary_ - T_))

        // Heat-transfer to the wall
      - fvm::Sp(htcw_->h()/Cp_, hs)
      + htcw_->h()*(hs/Cp_ + alpha_*(Tw_- T_))
    );
}


void thermoSingleLayer::solveEnergy()
{
    DebugInFunction << endl;

    updateSurfaceTemperatures();

    solve
    (
        fvm::ddt(deltaRho_, hs_)
      + fvm::div(phi_, hs_)
     ==
      - hsSp_
      - rhoSp_*hs_
      + q(hs_)
      + radiation_->Shs()
      + hsSpPyr_ // melt
    );

    correctThermoFields();

    // Evaluate viscosity from user-model
    viscosity_->correct(pPrimary_, T_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermoSingleLayer::thermoSingleLayer
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType,
    const bool readFields
)
:
    kinematicSingleLayer(modelType, mesh, g, regionType, false),
    thermo_(mesh.lookupObject<SLGThermo>("SLGThermo")),
    Cp_
    (
        IOobject
        (
            "Cp",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimTime/dimLength/dimTemperature, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),

    T_
    (
        IOobject
        (
            "Tf",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    Ts_
    (
        IOobject
        (
            "Tsf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T_,
        fvPatchFieldBase::zeroGradientType()
    ),
    Tw_
    (
        IOobject
        (
            "Twf",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        T_,
        fvPatchFieldBase::zeroGradientType()
    ),
    hs_
    (
        IOobject
        (
            "hf",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimMass, Zero),
        hsBoundaryTypes()
    ),
    qGasToFilm_
    (
     IOobject
     (
      "qGasToFilm",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
     zeroGradientFvPatchScalarField::typeName
    ),
    qFilmToWall_
    (
     IOobject
     (
      "qFilmToWall",
      time().timeName(),
      regionMesh(),
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     regionMesh(),
     dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
     zeroGradientFvPatchScalarField::typeName
    ),

    massAbs_
    (
        IOobject
        (
            "massAbsorption",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    energyAbs_
    (
        IOobject
        (
            "energyMassAbsorption",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    primaryEnergyTrans_
    (
        IOobject
        (
            "primaryEnergyTrans",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),

    charFrac_
    (
        IOobject
        (
            "charFrac",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    transmissivity_
    (
        IOobject
        (
            "filmTransmissivity",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("one", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    radTransmitted_
    (
        IOobject
        (
            "radTransmitted",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimArea/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    deltaWet_(coeffs_.get<scalar>("deltaWet")),
    hydrophilic_(coeffs_.get<bool>("hydrophilic")),
    hydrophilicDryScale_(0.0),
    hydrophilicWetScale_(0.0),

    hsSp_
    (
        IOobject
        (
            "hsSp",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimEnergy/dimArea/dimTime, Zero),
        this->mappedPushedFieldPatchTypes<scalar>()
    ),

    hsSpPrimary_
    (
        IOobject
        (
            hsSp_.name(), // Must have same name as hSp_ to enable mapping
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(hsSp_.dimensions(), Zero)
    ),
    impEnergyPrimary_
    (
        IOobject
        (
            "impEnergyPrimary", 
            time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        primaryMesh(),
        dimensionedScalar(dimEnergy, Zero)
    ),
    // Phase change (gasification) enthalpy
    hsSpGas_
    (
        IOobject
        (
            "hsSpGas",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(hsSp_.dimensions(), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    // Melt enthalpy
    hsSpPyr_
    (
        IOobject
        (
            "hsSpPyr",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(hsSp_.dimensions(), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    // Impingement enthalpy
    hsSpImp_
    (
        IOobject
        (
            "hsSpImp",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(hsSp_.dimensions(), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    hsSpDrip_
    (
        IOobject
        (
            "hsSpDrip",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(hsSp_.dimensions(), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    TPrimary_
    (
        IOobject
        (
            "T", // Same name as T on primary region to enable mapping
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimTemperature, Zero),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),
    YPrimary_(),
    // Additional heat flux
    massLossRate_
    (
        IOobject
        (
            "solidMLR",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime/dimArea, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    correctionZone_
    (
        IOobject
        (
            "correctionZone",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    burntZone_
    (
        IOobject
        (
            "burntZone",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            //IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    phiSolid_
    (
        IOobject
        (
            "mlrSolid", // Must match name defined in primary region (solver/include/infoFieldsOutput.H)
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimArea/dimTime, 0.0),
        this->mappedFieldAndInternalPatchTypes<scalar>()
    ),
    preHeatingCells_(coeffs_.getOrDefault<label>("preHeatingCells", 1)),
    burntMlrCrit_(coeffs_.getOrDefault<scalar>("burntMlrCrit", 1e-3)),
    viscosity_(filmViscosityModel::New(*this, coeffs(), mu_)),
    htcs_
    (
        heatTransferModel::New(*this, coeffs().subDict("upperSurfaceModels"))
    ),
    htcw_
    (
        heatTransferModel::New(*this, coeffs().subDict("lowerSurfaceModels"))
    ),
    phaseChange_(phaseChangeModel::New(*this, coeffs())),
    phaseChangeHocRatio_(coeffs().template
        getCheckOrDefault<scalar>("hocRatio",1,scalarMinMax::zero_one())),
    massAbsorption_(massAbsorptionModel::New(*this, coeffs())), // kvm
    radiation_(filmRadiationModel::New(*this, coeffs())),
    withTbounds_(limitType::CLAMP_NONE),
    Tbounds_(0, 5000)
{
    unsigned userLimits(limitType::CLAMP_NONE);

    if (coeffs().readIfPresent("Tmin", Tbounds_.min()))
    {
        userLimits |= limitType::CLAMP_MIN;
        Info<< "    limiting minimum temperature to " << Tbounds_.min() << nl;
    }

    if (coeffs().readIfPresent("Tmax", Tbounds_.max()))
    {
        userLimits |= limitType::CLAMP_MAX;
        Info<< "    limiting maximum temperature to " << Tbounds_.max() << nl;
    }
    withTbounds_ = limitType(userLimits);

    Info << "HOC conversion ratio:" << '\t' << phaseChangeHocRatio_ << endl;

    if (thermo_.hasMultiComponentCarrier())
    {
        YPrimary_.setSize(thermo_.carrier().species().size());

        forAll(thermo_.carrier().species(), i)
        {
            YPrimary_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        thermo_.carrier().species()[i],
                        time().timeName(),
                        regionMesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    regionMesh(),
                    dimensionedScalar(dimless, Zero),
                    pSp_.boundaryField().types()
                )
            );
        }
    }

    if (hydrophilic_)
    {
        coeffs_.readEntry("hydrophilicDryScale", hydrophilicDryScale_);
        coeffs_.readEntry("hydrophilicWetScale", hydrophilicWetScale_);
    }

    if (readFields)
    {
        transferPrimaryRegionThermoFields();

        correctAlpha();

        correctThermoFields();

        // Update derived fields
        hs_ == filmThermo_->hs(T_,true); // Alex

        deltaRho_ == delta_*rho_;

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

        // Evaluate viscosity from user-model
        viscosity_->correct(pPrimary_, T_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermoSingleLayer::~thermoSingleLayer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermoSingleLayer::addSources
(
    const label patchi,
    const label facei,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{
    kinematicSingleLayer::addSources
    (
        patchi,
        facei,
        massSource,
        momentumSource,
        pressureSource,
        energySource
    );

    DebugInfo
        << "    energy   = " << energySource << nl << nl;

    impEnergyPrimary_.boundaryFieldRef()[patchi][facei] -= energySource; // Alex, relax impingement sources

    //hsSpPrimary_.boundaryFieldRef()[patchi][facei] -= energySource;
}


void thermoSingleLayer::preEvolveRegion()
{
    DebugInFunction << endl;

    kinematicSingleLayer::preEvolveRegion();

    //- At this stage hsSp_ has only picked up impingement information 
    //- from transferPrimaryRegionSourceFields()
    //- Later, dripping and other stuff will be lumped in via updateSubmodels()
    //- So, the idea here is to cache the impingement information
    hsSpImp_ = hsSp_;
    hsSpImp_.correctBoundaryConditions();

    // Update phase change
    primaryMassTrans_ == dimensionedScalar("zero", dimMass, 0.0);
    primaryEnergyTrans_ == dimensionedScalar("zero", dimEnergy, 0.0);
    // Update mass absorption
    massAbs_ == dimensionedScalar("zero", dimMass, 0.0); // kvm
    energyAbs_ == dimensionedScalar("zero", dimEnergy, 0.0); // kvm
}


void thermoSingleLayer::evolveRegion()
{
    DebugInFunction << endl;    // Update film coverage indicator
    correctAlpha();

    // Solve continuity for deltaRho_
    solveContinuity();

    // Update sub-models to provide updated source contributions
    updateSubmodels();

    // Solve continuity for deltaRho_
    solveContinuity();

    for (int oCorr=1; oCorr<=nOuterCorr_; oCorr++)
    {
        // Explicit pressure source contribution
        tmp<volScalarField> tpu(this->pu());

        // Implicit pressure source coefficient
        tmp<volScalarField> tpp(this->pp());

        // Solve for momentum for U_
        tmp<fvVectorMatrix> tUEqn = solveMomentum(tpu(), tpp());
        fvVectorMatrix& UEqn = tUEqn.ref();

        // Solve energy for hs_ - also updates thermo
        solveEnergy();

        // Film thickness correction loop
        for (int corr=1; corr<=nCorr_; corr++)
        {
            // Solve thickness for delta_
            solveThickness(tpu(), tpp(), UEqn);
        }

        T_ == invertEnthalpy(); // Alex

    }


    // Update deltaRho_ with new delta_
    deltaRho_ == delta_*rho_;

    // Update film wall and surface velocities
    updateSurfaceVelocities(); // kvm

    // Update film wall and surface temperatures
    // updateSurfaceTemperatures();

    // Reset source terms for next time integration
    resetPrimaryRegionSourceTerms();
}


const volScalarField& thermoSingleLayer::Cp() const
{
    return Cp_;
}


const volScalarField& thermoSingleLayer::kappa() const
{
    return kappa_;
}


const volScalarField& thermoSingleLayer::T() const
{
    return T_;
}


const volScalarField& thermoSingleLayer::Ts() const
{
    return Ts_;
}


const volScalarField& thermoSingleLayer::Tw() const
{
    return Tw_;
}


const volScalarField& thermoSingleLayer::hs() const
{
    return hs_;
}


tmp<volScalarField> thermoSingleLayer::massAbs() const // kvm
{
    return massAbs_; // kvm
}


tmp<volScalarField> thermoSingleLayer::primaryMassTrans() const
{
    return primaryMassTrans_;
}


void thermoSingleLayer::info()
{
    kinematicSingleLayer::info();

    const scalarField& Tinternal = T_;

    Info<< indent << "min/mean/max(T)    = "
        << gMin(Tinternal) << ", "
        << gAverage(Tinternal) << ", "
        << gMax(Tinternal) << nl;

    phaseChange_->info(Info);
    massAbsorption_->info(Info); // kvm
}


tmp<volScalarField::Internal> thermoSingleLayer::Srho() const
{
    tmp<volScalarField::Internal> tSrho
    (
        new volScalarField::Internal
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
        )
    );

    scalarField& Srho = tSrho.ref();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs(), i)
    {
        const label filmPatchi = intCoupledPatchIDs()[i];

        scalarField patchMass =
            primaryMassTrans_.boundaryField()[filmPatchi];

        // hoc_solid/hoc_gas_species
        patchMass *= phaseChangeHocRatio_; 

        toPrimary(filmPatchi, patchMass);

        const label primaryPatchi = primaryPatchIDs()[i];
        const labelUList& cells =
            primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

        forAll(patchMass, j)
        {
            Srho[cells[j]] += patchMass[j]/(V[cells[j]]*dt);
        }
    }

    return tSrho;
}


tmp<volScalarField::Internal> thermoSingleLayer::Srho
(
    const label i
) const
{
    const label vapId = thermo_.carrierId(filmThermo_->name());

    tmp<volScalarField::Internal> tSrho
    (
        new volScalarField::Internal
        (
            IOobject
            (
                typeName + ":Srho(" + Foam::name(i) + ")",
                time_.timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            primaryMesh(),
            dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
        )
    );

    if (vapId == i)
    {
        scalarField& Srho = tSrho.ref();
        const scalarField& V = primaryMesh().V();
        const scalar dt = time().deltaTValue();

        forAll(intCoupledPatchIDs_, i)
        {
            const label filmPatchi = intCoupledPatchIDs_[i];

            scalarField patchMass =
                primaryMassTrans_.boundaryField()[filmPatchi];

            patchMass *= phaseChangeHocRatio_; 

            toPrimary(filmPatchi, patchMass);

            const label primaryPatchi = primaryPatchIDs()[i];
            const labelUList& cells =
                primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

            forAll(patchMass, j)
            {
                Srho[cells[j]] += patchMass[j]/(V[cells[j]]*dt);
            }
        }
    }

    return tSrho;
}


tmp<volScalarField::Internal> thermoSingleLayer::Sh() const
{
    tmp<volScalarField::Internal> tSh
    (
        new volScalarField::Internal
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
        )
    );
/*
    phase change energy fed back into the film...

    scalarField& Sh = tSh.ref();
    const scalarField& V = primaryMesh().V();
    const scalar dt = time_.deltaTValue();

    forAll(intCoupledPatchIDs_, i)
    {
        const label filmPatchi = intCoupledPatchIDs_[i];

        scalarField patchEnergy =
            primaryEnergyTrans_.boundaryField()[filmPatchi];

        toPrimary(filmPatchi, patchEnergy);

        const label primaryPatchi = primaryPatchIDs()[i];
        const labelUList& cells =
            primaryMesh().boundaryMesh()[primaryPatchi].faceCells();

        forAll(patchEnergy, j)
        {
            Sh[cells[j]] += patchEnergy[j]/(V[cells[j]]*dt);
        }
    }
*/
    return tSh;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace regionModels
} // End namespace surfaceFilmModels

// ************************************************************************* //
