/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "inertPlateOneSideHeatingFvPatchScalarField.H"
#include "LESModel.H"
#include "basicThermo.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
#include "turbulentFluidThermoModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inertPlateOneSideHeatingFvPatchScalarField::
inertPlateOneSideHeatingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    density_(7800.0),
    Cp_(500),
    thickness_(0.002),
    includeRadiation_(false),
    backRadiationLoss_(false),
    heatLoss_(10.0)
{
    //checkType();
    //read();
}


inertPlateOneSideHeatingFvPatchScalarField::
inertPlateOneSideHeatingFvPatchScalarField
(
    const inertPlateOneSideHeatingFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    density_(ptf.density_),
    Cp_(ptf.Cp_),
    thickness_(ptf.thickness_),
    includeRadiation_(ptf.includeRadiation_),
    backRadiationLoss_(ptf.backRadiationLoss_),
    heatLoss_(ptf.heatLoss_)
{}


inertPlateOneSideHeatingFvPatchScalarField::
inertPlateOneSideHeatingFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    density_(dict.lookupOrDefault<scalar>("density", 7800.0)),
    Cp_(dict.lookupOrDefault<scalar>("Cp", 500)),
    thickness_(dict.lookupOrDefault<scalar>("thickness", 0.002)),
    includeRadiation_(dict.lookupOrDefault<bool>("includeRadiation", false)),
    backRadiationLoss_(dict.lookupOrDefault<bool>("backSideRadiationLoss", false)),
    heatLoss_(dict.lookupOrDefault<scalar>("heatLoss", 10.0))
{
//    checkType();
    //read();
}


inertPlateOneSideHeatingFvPatchScalarField::
inertPlateOneSideHeatingFvPatchScalarField
(
    const inertPlateOneSideHeatingFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    density_(tppsf.density_),
    Cp_(tppsf.Cp_),
    thickness_(tppsf.thickness_),
    includeRadiation_(tppsf.includeRadiation_),
    backRadiationLoss_(tppsf.backRadiationLoss_),
    heatLoss_(tppsf.heatLoss_)
{
//    checkType();
}


inertPlateOneSideHeatingFvPatchScalarField::
inertPlateOneSideHeatingFvPatchScalarField
(
    const inertPlateOneSideHeatingFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    density_(tppsf.density_),
    Cp_(tppsf.Cp_),
    thickness_(tppsf.thickness_),
    includeRadiation_(tppsf.includeRadiation_),
    backRadiationLoss_(tppsf.backRadiationLoss_),
    heatLoss_(tppsf.heatLoss_)
{
//    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void inertPlateOneSideHeatingFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    const scalarField& qrW = patch().lookupPatchField<volScalarField, scalar>("qr");
    //const scalarField& qcW = patch().lookupPatchField<volScalarField, scalar>("QcWallFunction");
    const scalarField& qcWT = patch().lookupPatchField<surfaceScalarField, scalar>("convectiveHeatFlux_T");

    volVectorField sprayFlux = db().lookupObject<volVectorField>("sprayFluxNet");
    const labelUList& faceCells = patch().faceCells();

    scalarField& Twall = *this;

    forAll(Twall,facei)
    {
        scalar fluxZ = max(0, sprayFlux[faceCells[facei]].z());
        scalar qWater(0);
        if(Twall[facei] > 373)
        {
            qWater = fluxZ*2.6e6;
        }
        else
        {
            qWater = fluxZ*(Twall[facei]-298)*4.2e3;
        }

        scalar qWall(qcWT[facei]);
        if(includeRadiation_)
        {
            qWall = qWall + qrW[facei];
        }
        if(backRadiationLoss_)
        {
	    qWall = qWall - constant::physicoChemical::sigma.value()*(pow(Twall[facei],4) - pow(300.0,4));
        }

        scalar dT = db().time().deltaT().value()
                   *(qWall - heatLoss_*(Twall[facei]-300.0) - qWater)
                   /(density_*Cp_*thickness_);

        Twall[facei] =min(1500,max(293, Twall[facei]+dT));
    }

}

void inertPlateOneSideHeatingFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    //writeLocalEntries(os);
    os.writeKeyword("density") << density_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cp") << Cp_ << token::END_STATEMENT << nl;
    os.writeKeyword("thickness") << thickness_ << token::END_STATEMENT << nl;
    os.writeKeyword("heatLoss") << heatLoss_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    inertPlateOneSideHeatingFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
