/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2023 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "alphatPyrolysisMeshCorrectedFvPatchScalarField.H"
#include "LESModel.H"
#include "basicThermo.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatPyrolysisMeshCorrectedFvPatchScalarField::
alphatPyrolysisMeshCorrectedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    phiName_("phi"),
    fuelConversionRatio_(2.5), // ratio of hoc_gasfuel/hoc_pyr
    m60_(0.01),  // normalizing factor to set blowing factor to 0.6 if mlr is 10 g/m2/s, reported by J. deRis
    Prt_(1.0),  // turbulent prandlt number
    C1_(175.0), // 1st (linear) coefficient of polynomial for imposing grid-size independence
    C2_(0), // 2nd (parabolic) coefficient of polynomial for imposing grid-size independence
    LR_(0.0015), // If mesh size less than this value, considered fully-resolved
    LB_(0.0085) // Length-scale of blowing effect
{

}


alphatPyrolysisMeshCorrectedFvPatchScalarField::
alphatPyrolysisMeshCorrectedFvPatchScalarField
(
    const alphatPyrolysisMeshCorrectedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    fuelConversionRatio_(ptf.fuelConversionRatio_),
    m60_(ptf.m60_),
    Prt_(ptf.Prt_),
    C1_(ptf.C1_),
    C2_(ptf.C2_),
    LR_(ptf.LR_),
    LB_(ptf.LB_)
{}


alphatPyrolysisMeshCorrectedFvPatchScalarField::
alphatPyrolysisMeshCorrectedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    phiName_(dict.lookupOrDefault<word>("phiName","phi")),
    fuelConversionRatio_(dict.lookupOrDefault<scalar>("fuelConversionRatio", 2.5)), 
    m60_(dict.lookupOrDefault<scalar>("m60", 0.01)),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 1.0)),
    C1_(dict.lookupOrDefault<scalar>("C1", 225.0)),
    C2_(dict.lookupOrDefault<scalar>("C2", 0)),
    LR_(dict.lookupOrDefault<scalar>("LR", 0.0015)),
    LB_(dict.lookupOrDefault<scalar>("LB",0.0085))
{
}


alphatPyrolysisMeshCorrectedFvPatchScalarField::
alphatPyrolysisMeshCorrectedFvPatchScalarField
(
    const alphatPyrolysisMeshCorrectedFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    phiName_(tppsf.phiName_),
    fuelConversionRatio_(tppsf.fuelConversionRatio_),
    m60_(tppsf.m60_),
    Prt_(tppsf.Prt_),
    C1_(tppsf.C1_),
    C2_(tppsf.C2_),
    LR_(tppsf.LR_),
    LB_(tppsf.LB_)
{
}


alphatPyrolysisMeshCorrectedFvPatchScalarField::
alphatPyrolysisMeshCorrectedFvPatchScalarField
(
    const alphatPyrolysisMeshCorrectedFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    phiName_(tppsf.phiName_),
    fuelConversionRatio_(tppsf.fuelConversionRatio_),
    m60_(tppsf.m60_),
    Prt_(tppsf.Prt_),
    C1_(tppsf.C1_),
    C2_(tppsf.C2_),
    LR_(tppsf.LR_),
    LB_(tppsf.LB_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatPyrolysisMeshCorrectedFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    const label patchi = patch().index();

    const compressible::turbulenceModel& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const basicThermo& thermo = db().lookupObject<basicThermo>
    (
        "thermophysicalProperties"
    );

    const fvPatchScalarField& alphaw = turbModel.alpha()().boundaryField()[patchi];
    scalarField& alphatw = *this;

    const fvPatchScalarField& Tw = thermo.T().boundaryField()[patchi];
    const scalarField T(Tw.patchInternalField());

    const labelUList& faceCells = patch().faceCells();
    const scalarField& yw = turbModel.y()[patchi];

    // Assign gas flux to phiw based on name and field type
    // Usually phiName_ will correspond to surfaceScalarField "phi"
    // But named field may of type volScalarField, handle that here
    scalarField phiw(patch().size(),Zero);

    if (this->db().template foundObject<surfaceScalarField>(phiName_))
    {
        phiw = patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);
    }
    else if(this->db().template foundObject<volScalarField>(phiName_))
    {
        // Volume scalar field value is positive
        // Later definition of pyrolysisMassFlux assumes phiw is negative
        // as this is the case for the flux term from the pyrolysis region
        // In this code-block, the lookup of phiName_ should return positive 
        // value. 
        phiw = -patch().lookupPatchField<volScalarField, scalar>(phiName_);
    }

    // phi value from coupled pyrolysis region is negative
    scalarField pyrolysisMassFlux( - fuelConversionRatio_*phiw/patch().magSf());

    forAll(alphatw,facei)
    {
        /*
            1. Mesh size correction factor
        */

        // Near wall resolution deficit parameter
        scalar dDelta(max(0,yw[facei]*2.0-LR_));

        // Equation 2-4 (mesh correction) 
        // truncated to second order term
        scalar Mc(1 + C1_*dDelta + C2_*dDelta*dDelta);

        /*
            2. Pyrolysis (blowing effect) correction factor
        */

        // mRatio is a measure of how much blowing there is relative to
        // the value associated with a 60% reduction in convective heat flux
        // corresponds to  mdot``/mdot``_0 term shown in Equation 2-6
        scalar mRatio(pyrolysisMassFlux[facei]/m60_ + 0.001);

        // Equation 2-6 
        scalar BF(mRatio/(exp(mRatio)-1.0));

        // Equation 2-5 
        scalar Bc(BF + (1-BF)*exp(-dDelta/LB_));

        /*
            3. Evaluate correction alphat
        */
        
        // Equation 2-3 made explicit and with the additional correction terms
        alphatw[facei] = alphaw[facei]*(Bc*Mc - 1.0);

    }
}

void alphatPyrolysisMeshCorrectedFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("phiName") << phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fuelConversionRatio") << fuelConversionRatio_ << token::END_STATEMENT << nl;
    os.writeKeyword("m60") << m60_ << token::END_STATEMENT << nl;
    os.writeKeyword("Prt") << Prt_ << token::END_STATEMENT << nl;
    os.writeKeyword("C1") << C1_ << token::END_STATEMENT << nl;
    os.writeKeyword("C2") << C2_ << token::END_STATEMENT << nl;
    os.writeKeyword("LR") << LR_ << token::END_STATEMENT << nl;
    os.writeKeyword("LB") << LB_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatPyrolysisMeshCorrectedFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
