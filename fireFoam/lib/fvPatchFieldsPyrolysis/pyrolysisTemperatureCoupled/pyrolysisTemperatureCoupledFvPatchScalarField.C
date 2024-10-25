/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "pyrolysisTemperatureCoupledFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "mapDistribute.H"
#include "regionProperties.H"
#include "basicThermo.H"
#include "turbulenceModel.H"

#include "constants.H"

#include "radiationModel.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pyrolysisTemperatureCoupledFvPatchScalarField::
pyrolysisTemperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase
    (
        patch()
    ),
    pyrolysisRegionName_("pyrolysisProperties"),
    TnbrName_("T"),
    neighbourFieldRadiativeName_("qin"),
    neighbourFieldConvectiveName_("convectiveHeatFlux_T"),
    KName_("K")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


pyrolysisTemperatureCoupledFvPatchScalarField::
pyrolysisTemperatureCoupledFvPatchScalarField
(
    const pyrolysisTemperatureCoupledFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    pyrolysisRegionName_(psf.pyrolysisRegionName_),
    TnbrName_(psf.TnbrName_),
    neighbourFieldRadiativeName_(psf.neighbourFieldRadiativeName_),
    neighbourFieldConvectiveName_(psf.neighbourFieldConvectiveName_),
    KName_(psf.KName_)
{}


pyrolysisTemperatureCoupledFvPatchScalarField::
pyrolysisTemperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    pyrolysisRegionName_
    (
        dict.lookupOrDefault<word>("pyrolysisRegion", "pyrolysisProperties")
    ),
    TnbrName_(dict.lookupOrDefault<word>("Tnbr","T")),
    neighbourFieldRadiativeName_(dict.lookupOrDefault<word>("neighbourFieldRadiativeName","qin")),
    neighbourFieldConvectiveName_(dict.lookupOrDefault<word>("neighbourFieldConvectiveName","convectiveHeatFlux_T")),
    KName_(dict.lookupOrDefault<word>("K","K"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "pyrolysisTemperatureCoupledFvPatchScalarField::"
            "pyrolysisTemperatureCoupledFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


pyrolysisTemperatureCoupledFvPatchScalarField::
pyrolysisTemperatureCoupledFvPatchScalarField
(
    const pyrolysisTemperatureCoupledFvPatchScalarField&
        psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    pyrolysisRegionName_(psf.pyrolysisRegionName_),
    TnbrName_(psf.TnbrName_),
    neighbourFieldRadiativeName_(psf.neighbourFieldRadiativeName_),
    neighbourFieldConvectiveName_(psf.neighbourFieldConvectiveName_),
    KName_(psf.KName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pyrolysisTemperatureCoupledFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void pyrolysisTemperatureCoupledFvPatchScalarField::rmap
(
    const fvPatchScalarField& psf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(psf, addr);
}



void pyrolysisTemperatureCoupledFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    const label patchI = patch().index();
    const label nbrPatchI = mpp.samplePolyPatch().index();
    const polyMesh& mesh = patch().boundaryMesh().mesh();
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[nbrPatchI];

    scalarField
        nbrField =
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
        );


    mpp.distribute(nbrField);

    const scalarField K(this->kappa(*this));

    scalarList radiField(nbrPatch.size(), 0.0);
    scalarList convField(nbrPatch.size(), 0.0);

    const scalar sigma = constant::physicoChemical::sigma.value();

    // In solid
    if(mesh.name() == pyrolysisRegionName_) 
    {
        const fvMesh& mesh = patch().boundaryMesh().mesh();
         if (! (mesh.foundObject<radiation::radiationModel>("radiationProperties")))
         {
             FatalErrorIn
             (
                 "pyrolysisTemperatureCoupledFvPatchScalarField::"
                 "pyrolysisTemperatureCoupledFvPatchScalarField\n"
                 "(\n"
                 "    const fvPatch& p,\n"
                 "    const DimensionedField<scalar, volMesh>& iF,\n"
                 "    const dictionary& dict\n"
                 ")\n"
             )   << "\n    radiationProperties file not found in pyrolysis region\n"
                 << exit(FatalError);
         }
         const radiation::radiationModel& radiation =
             mesh.lookupObject<radiation::radiationModel>
             (
                 "radiationProperties"
             );

         scalarField temissivity
         (
             radiation.absorptionEmission().e()().boundaryField()
             [
                 patch().index()
             ]
         );

         scalarField tabsorptivity
         (
             radiation.absorptionEmission().a()().boundaryField()
             [
                 patch().index()
             ]
         );

        radiField =
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                neighbourFieldRadiativeName_
            );

	    mpp.distribute(radiField);
        convField =	
            nbrPatch.lookupPatchField<surfaceScalarField, scalar>
            (
                neighbourFieldConvectiveName_
            );

	    mpp.distribute(convField);

        forAll(*this, i)
        {
            scalar qConv = -convField[i];
            scalar qRad  = -tabsorptivity[i]*radiField[i] + temissivity[i]*sigma*pow4(operator[](i));
            scalar qNet  = qConv + qRad;

            this->refValue()[i] = operator[](i); 
            this->refGrad()[i] = -qNet/K[i];

            this->valueFraction()[i] = 0.0; // Neumann 
        }
    }
    else if (mesh.name() == "region0") // gas region 
    {
        this->refValue() = nbrField;
        this->refGrad() = 0.0;   
        this->valueFraction() = 1.0; // Dirichlet
    }
    else
    {
        FatalErrorIn
        (
            "const pyrolysisTemperatureCoupledFvPatchScalarField::"
            "updateCoeffs() "
        )
        << "Ill-defined regionName in BC specification, got:  " << mesh.name() 
        << abort(FatalError);
    }


    mixedFvPatchScalarField::updateCoeffs();
}


void pyrolysisTemperatureCoupledFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>
    (
        "pyrolysisRegion",
        "pyrolysisProperties",
        pyrolysisRegionName_
    );
    os.writeKeyword("Tnbr")<< TnbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldRadiativeName")<<
        neighbourFieldRadiativeName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldConvectiveName")<< 
        neighbourFieldConvectiveName_ << token::END_STATEMENT << nl;
    os.writeKeyword("K")<< 
        KName_ << token::END_STATEMENT << nl;
    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    pyrolysisTemperatureCoupledFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
