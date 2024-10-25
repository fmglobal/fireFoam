/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "pyrolysisVelocityCoupledFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "mapDistribute.H"
#include "basicThermo.H"
#include "surfaceFields.H"

#include "singleStepReactingMixture.H"
#include "thermoPhysicsTypes.H"

#include "reactingMixture.H"
#include "constIsoSolidTransport.H"
#include "hConstThermo.H"
#include "rhoConst.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pyrolysisVelocityCoupledFvPatchVectorField::
pyrolysisVelocityCoupledFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    nbrPhiName_("none"),
    phiName_("phi"),
    rhoName_("rho"),
    virginName_("v"),
    charName_("char"),
    hocVirgin_(0.0),
    hocChar_(0.0),
    charring_(false),
    complexFuel_(false)
{}


Foam::pyrolysisVelocityCoupledFvPatchVectorField::
pyrolysisVelocityCoupledFvPatchVectorField
(
    const pyrolysisVelocityCoupledFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    nbrPhiName_(ptf.nbrPhiName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    virginName_(ptf.virginName_),
    charName_(ptf.charName_),
    hocVirgin_(ptf.hocVirgin_),
    hocChar_(ptf.hocChar_),
    charring_(ptf.charring_),
    complexFuel_(ptf.complexFuel_)
{}


Foam::pyrolysisVelocityCoupledFvPatchVectorField::
pyrolysisVelocityCoupledFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    nbrPhiName_(dict.getOrDefault<word>("nbrPhi", "phi")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    virginName_(dict.getOrDefault<word>("virgin","v")),
    charName_(dict.getOrDefault<word>("char","char")),
    hocVirgin_(dict.get<scalar>("hocVirgin")),
    hocChar_(dict.getOrDefault<scalar>("hocChar",32.8e6)),
    charring_(dict.get<bool>("charring")),
    complexFuel_(dict.getOrDefault<bool>("complexFuel", false))
{}


Foam::pyrolysisVelocityCoupledFvPatchVectorField::
pyrolysisVelocityCoupledFvPatchVectorField
(
    const pyrolysisVelocityCoupledFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    nbrPhiName_(ptf.nbrPhiName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    virginName_(ptf.virginName_),
    charName_(ptf.charName_),
    hocVirgin_(ptf.hocVirgin_),
    hocChar_(ptf.hocChar_),
    charring_(ptf.charring_),
    complexFuel_(ptf.complexFuel_)
{}


Foam::pyrolysisVelocityCoupledFvPatchVectorField::
pyrolysisVelocityCoupledFvPatchVectorField
(
    const pyrolysisVelocityCoupledFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    nbrPhiName_(ptf.nbrPhiName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    virginName_(ptf.virginName_),
    charName_(ptf.charName_),
    hocVirgin_(ptf.hocVirgin_),
    hocChar_(ptf.hocChar_),
    charring_(ptf.charring_),
    complexFuel_(ptf.complexFuel_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pyrolysisVelocityCoupledFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        patch().patch()
    );
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch = refCast<const fvMesh>
    (
        nbrMesh
    ).boundary()[mpp.samplePolyPatch().index()];

    scalarList phi =
        nbrPatch.lookupPatchField<surfaceScalarField, scalar>(nbrPhiName_);

    // get heat of combustion of the gasious fuel
    const basicThermo& thermo =
        db().lookupObject<basicThermo>("thermophysicalProperties");

    const singleStepReactingMixture<gasHThermoPhysics>& singleMixture
    (
        dynamic_cast<const singleStepReactingMixture<gasHThermoPhysics>&>
        (thermo)
    );

    // heat of combustion [J/kg]
    scalar qFuel(singleMixture.qFuel().value());

    scalar hocPyr(hocVirgin_);

    if (charring_)
    {
        const basicThermo& solidThermo =
         mpp.sampleMesh().lookupObject<basicThermo>("thermophysicalProperties");

        // access density of char and v
        const reactingMixture<constIsoSolidTransport<species::thermo<hConstThermo<rhoConst<specie> >, sensibleEnthalpy > > >& solidMixture
        (
            dynamic_cast<const reactingMixture<constIsoSolidTransport<species::thermo<hConstThermo<rhoConst<specie> >, sensibleEnthalpy > > >&>
            (solidThermo)
        );

        const label charIndex = solidMixture.species()[charName_];
        const label vIndex = solidMixture.species()[virginName_];

        const scalar rhoChar = solidMixture.rho(charIndex,1,298.0);  //plug in any P and T for rhoConst
        const scalar rhoV = solidMixture.rho(vIndex,1,298.0); 

        // Heat of combustion of gaseous pyrolysate
        // The input value hocVirgin_ is the measured HOC including all combustible content 
        // We assume that charring_ commodities leave an unburnt char layer, the HOC associated with
        // the char must therefore be deducted from the measured value to get the effective value
        // The numerator is the chemical energy liberated from the solid (excluding char contribution)
        // The denominator is the mass of gas produced per mass of virgin solid.
        scalar hocPyr = (hocVirgin_ * rhoV - hocChar_ * rhoChar) / (rhoV - rhoChar);
    }

    scalarField hocField(nbrPatch.size(),hocPyr);

    // For CUP model only, switches between cardboard liner and CUP HOC
    // CUP model responsible for setting the hocComplex value, which can
    // vary dynamically
    if (complexFuel_)
    {
        const scalarField& isLinerGone =
            nbrPatch.lookupPatchField<volScalarField, scalar>("swccGone");
        const scalarField& hocComplex =
            nbrPatch.lookupPatchField<volScalarField, scalar>("pyrolHOC");
        const scalar& hocPyrCC = hocPyr;

        hocField  = 0.5*(1. + isLinerGone) * hocComplex + 0.5*(1. - isLinerGone) * hocPyr;
    }

    // convert to equivalent gaseous fuel
    phi = phi * hocField / qFuel;

    mpp.distribute(phi);

    const surfaceScalarField& phiName =
        db().lookupObject<surfaceScalarField>(phiName_);

    scalarField U(-phi/patch().magSf());

    vectorField n(patch().nf());

    if (phiName.dimensions() == dimVolume/dimTime)
    {
        // volumetric flow-rate
        operator==(n*U);
    }
    else if (phiName.dimensions() == dimMass/dimTime)
    {
        const auto& rhop =
            patch().lookupPatchField<volScalarField>(rhoName_);

        // mass flow-rate
        operator==(n*U/rhop);

        if (debug)
        {
            scalar phi(gSum(rhop*(*this) & patch().Sf()));
            Info<< patch().boundaryMesh().mesh().name() << ':'
                << patch().name() << ':'
                << this->internalField().name() << " <- "
                << nbrMesh.name() << ':'
                << nbrPatch.name() << ':'
                << this->internalField().name() << " :"
                << " mass flux[Kg/s]:" << -phi
                << endl;
        }
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of " << phiName_ << " are incorrect" << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << nl << exit(FatalError);
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::pyrolysisVelocityCoupledFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("virgin", "v", virginName_);
    os.writeEntryIfDifferent<word>("char", "char", charName_);
    os.writeKeyword("nbrPhi") << nbrPhiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("hocChar") << hocChar_ << token::END_STATEMENT << nl;
    os.writeKeyword("hocVirgin") << hocVirgin_ << token::END_STATEMENT << nl;
    os.writeKeyword("charring") << charring_ << token::END_STATEMENT << nl;
    os.writeKeyword("complexFuel") << complexFuel_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       pyrolysisVelocityCoupledFvPatchVectorField
   );
}


// ************************************************************************* //
