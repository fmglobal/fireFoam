/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenCFD Ltd.
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

#include "filmPyrolysisTemperatureCoupledFFFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "mapDistribute.H"
#include "regionProperties.H"
#include "thermoSingleLayer.H"
#include "pyrolysisModel.H"

#include "constants.H"

#include "radiationModel.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const filmPyrolysisTemperatureCoupledFFFvPatchScalarField::filmModelType&
filmPyrolysisTemperatureCoupledFFFvPatchScalarField::
filmModel() const
{
    HashTable<const filmModelType*> models
        = db().time().lookupClass<filmModelType>();

    forAllConstIter(HashTable<const filmModelType*>, models, iter)
    {
        if (iter()->regionMesh().name() == filmRegionName_)
        {
            return *iter();
        }
    }


    FatalErrorIn
    (
        "const filmPyrolysisTemperatureCoupledFFFvPatchScalarField::"
        "filmModelType& "
        "filmPyrolysisTemperatureCoupledFFFvPatchScalarField::"
        "filmModel() const"
    )
        << "Unable to locate film region " << filmRegionName_
        << abort(FatalError);

    return **models.begin();
}


const filmPyrolysisTemperatureCoupledFFFvPatchScalarField::
pyrolysisModelType&
filmPyrolysisTemperatureCoupledFFFvPatchScalarField::
pyrModel() const
{
    HashTable<const pyrolysisModelType*> models =
        db().time().lookupClass<pyrolysisModelType>();

    forAllConstIter(HashTable<const pyrolysisModelType*>, models, iter)
    {
        if (iter()->regionMesh().name() == pyrolysisRegionName_)
        {
            return *iter();
        }
    }

    FatalErrorIn
    (
        "const filmPyrolysisTemperatureCoupledFFFvPatchScalarField::"
        "pyrolysisModelType& "
        "filmPyrolysisTemperatureCoupledFFFvPatchScalarField::"
        "pyrModel() const"
    )
        << "Unable to locate pyrolysis region " << pyrolysisRegionName_
        << abort(FatalError);

    return **models.begin();
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

filmPyrolysisTemperatureCoupledFFFvPatchScalarField::
filmPyrolysisTemperatureCoupledFFFvPatchScalarField
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
    filmRegionName_("surfaceFilmProperties"),
    pyrolysisRegionName_("pyrolysisProperties"),
    TnbrName_("undefined-Tnbr"),
    neighbourFieldRadiativeName_("undefined-neigbourFieldRadiativeName"),
    neighbourFieldConvectiveName_("undefined-neigbourFieldConvectiveName"),
    fieldRadiativeName_("undefined-fieldRadiativeName"),
    fieldConvectiveName_("undefined-fieldConvectiveName"),
    KName_("undefined-K"),
    qExtra_(0),
    filmDeltaDry_(0),
    filmDeltaWet_(0),
    qRadiativeTransfer_(p.size(),0),
    qConvectiveTransfer_(p.size(),0),
    qExtraTransfer_(p.size(),0),
    wetnessFactorTransfer_(p.size(),0),
    filmAlphaTransfer_(p.size(),0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


filmPyrolysisTemperatureCoupledFFFvPatchScalarField::
filmPyrolysisTemperatureCoupledFFFvPatchScalarField
(
    const filmPyrolysisTemperatureCoupledFFFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    filmRegionName_(psf.filmRegionName_),
    pyrolysisRegionName_(psf.pyrolysisRegionName_),
    TnbrName_(psf.TnbrName_),
    neighbourFieldRadiativeName_(psf.neighbourFieldRadiativeName_),
    neighbourFieldConvectiveName_(psf.neighbourFieldConvectiveName_),
    fieldRadiativeName_(psf.fieldRadiativeName_),
    fieldConvectiveName_(psf.fieldConvectiveName_),
    KName_(psf.KName_),
    qExtra_(psf.qExtra_),
    filmDeltaDry_(psf.filmDeltaDry_),
    filmDeltaWet_(psf.filmDeltaWet_),
    qRadiativeTransfer_(psf.qRadiativeTransfer_),
    qConvectiveTransfer_(psf.qConvectiveTransfer_),
    qExtraTransfer_(psf.qConvectiveTransfer_),
    wetnessFactorTransfer_(psf.wetnessFactorTransfer_),
    filmAlphaTransfer_(psf.filmAlphaTransfer_)
{}


filmPyrolysisTemperatureCoupledFFFvPatchScalarField::
filmPyrolysisTemperatureCoupledFFFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    filmRegionName_
    (
        dict.lookupOrDefault<word>("filmRegion", "surfaceFilmProperties")
    ),
    pyrolysisRegionName_
    (
        dict.lookupOrDefault<word>("pyrolysisRegion", "pyrolysisProperties")
    ),
    TnbrName_(dict.lookup("Tnbr")),
    neighbourFieldRadiativeName_(dict.lookup("neighbourFieldRadiativeName")),
    neighbourFieldConvectiveName_(dict.lookup("neighbourFieldConvectiveName")),
    fieldRadiativeName_(dict.lookup("fieldRadiativeName")),
    fieldConvectiveName_(dict.lookup("fieldConvectiveName")),
    KName_(dict.lookup("K")),
    qExtra_(dict.lookupOrDefault<scalar>("qExtra", 0.0)),
    filmDeltaDry_(dict.lookupOrDefault<scalar>("filmDeltaDry",0)),
    filmDeltaWet_(dict.lookupOrDefault<scalar>("filmDeltaWet",0.0002)),
    qRadiativeTransfer_(p.size(),0),
    qConvectiveTransfer_(p.size(),0),
    qExtraTransfer_(p.size(),0),
    wetnessFactorTransfer_(p.size(),0),
    filmAlphaTransfer_(p.size(),0)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "filmPyrolysisTemperatureCoupledFFFvPatchScalarField::"
            "filmPyrolysisTemperatureCoupledFFFvPatchScalarField\n"
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


filmPyrolysisTemperatureCoupledFFFvPatchScalarField::
filmPyrolysisTemperatureCoupledFFFvPatchScalarField
(
    const filmPyrolysisTemperatureCoupledFFFvPatchScalarField&
        psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    filmRegionName_(psf.filmRegionName_),
    pyrolysisRegionName_(psf.pyrolysisRegionName_),
    TnbrName_(psf.TnbrName_),
    neighbourFieldRadiativeName_(psf.neighbourFieldRadiativeName_),
    neighbourFieldConvectiveName_(psf.neighbourFieldConvectiveName_),
    fieldRadiativeName_(psf.fieldRadiativeName_),
    fieldConvectiveName_(psf.fieldConvectiveName_),
    KName_(psf.KName_),
    qExtra_(psf.qExtra_),
    filmDeltaDry_(psf.filmDeltaDry_),
    filmDeltaWet_(psf.filmDeltaWet_),
    qRadiativeTransfer_(psf.qRadiativeTransfer_),
    qConvectiveTransfer_(psf.qConvectiveTransfer_),
    qExtraTransfer_(psf.qExtraTransfer_),
    wetnessFactorTransfer_(psf.wetnessFactorTransfer_),
    filmAlphaTransfer_(psf.filmAlphaTransfer_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void filmPyrolysisTemperatureCoupledFFFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void filmPyrolysisTemperatureCoupledFFFvPatchScalarField::rmap
(
    const fvPatchScalarField& psf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(psf, addr);
}


void filmPyrolysisTemperatureCoupledFFFvPatchScalarField::
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

    // This is the patch field from the 'other' patch, nbrField is the patch temperature
    //  - E.g., if calling from gas region this is the pyrolysis boundary value
    scalarField
        nbrField =
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
        );

    mpp.distribute(nbrField);

    // Thermal conductivity, used to match heat conduction into solid
    // from surface with net heat flux to solid
    const scalarField K(this->kappa(*this));

    scalarList nbrTotalFlux (nbrPatch.size(),0.0);
    mpp.distribute(nbrTotalFlux);

    scalarList filmRadTransmitted(nbrPatch.size(),0.0);
    scalarList filmTransmissivity(nbrPatch.size(),0.0);

    scalarList Tfilm(nbrPatch.size(), 0.0);
    scalarList alpha(nbrPatch.size(), 0.0);
    scalarList filmConv(nbrPatch.size(), 0.0);
    scalarList qExtraZone(nbrPatch.size(), 0.0);
    
    // Only used for CUP model
    scalarList filmDelta(nbrPatch.size(),0.0);
    
    scalarList convField(nbrPatch.size(), 0.0);

    const pyrolysisModelType& pyrolysis = pyrModel();

    const filmModelType& film = filmModel();

    if(mesh.name() == pyrolysisRegionName_) 
    {

        // Properties aquired from film model
        filmConv =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "qFilmToWall",
                patchI,
                true
            );

        filmRadTransmitted = 
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "radTransmitted",
                patchI,
                true
            );

        filmTransmissivity = 
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "filmTransmissivity",
                patchI,
                true
            );

        Tfilm =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "Tf",
                patchI,
                true
            );
        alpha =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "alpha",
                patchI,
                true
            );
        // Only used for CUP model
        filmDelta =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "deltaf",
                patchI,
                true
            );

        // Lateral flame spread term 
        // apply additional heat flux at pyrolysis front
        bool qExtraZoneExists = false;
        if(film.regionMesh().foundObject<volScalarField>("correctionZone"))
        {
            qExtraZone =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "correctionZone",
                patchI,
                true
            );
            Info<<"Add extra heat flux at MLR front: "<<qExtra_<<endl;
            qExtraZoneExists = true;
        }

       
        // Get radiation properties 
        if (! (mesh.foundObject<radiation::radiationModel>("radiationProperties")))
        {
            FatalErrorIn
            (
                "filmPyrolysisTemperatureCoupledFFFvPatchScalarField::"
                "filmPyrolysisTemperatureCoupledFFFvPatchScalarField\n"
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
                //nbrPatch.index()
                patch().index()
            ]
        );

        scalarField tabsorptivity
        (
            radiation.absorptionEmission().a()().boundaryField()
            [
                //nbrPatch.index()
                patch().index()
            ]
        );

        // Convective heat flux from the fluid region
        convField =	
            nbrPatch.lookupPatchField<surfaceScalarField, scalar>
            (
                neighbourFieldConvectiveName_
            );
 
	    mpp.distribute(convField);

        // Estimate wetness of the film (1: wet , 0: dry)
        // AK - Needed in complexFuelPyrolysis only (CUP)
        scalarField wetnessFraction((filmDelta - filmDeltaDry_)/max(SMALL, (filmDeltaWet_ - filmDeltaDry_)));
        wetnessFraction = min(scalar(1.), max(scalar(0.), wetnessFraction));

        forAll(*this, i)
        {
            scalar qConvWet = -filmConv[i];
            scalar qConvDry = -convField[i]; 
            // Blending function based on wetness fraction
            scalar qConv = alpha[i]*qConvWet + (1.0-alpha[i])*qConvDry;

            const scalar sigma = constant::physicoChemical::sigma.value();

            // - filmRadTransmitted accounts for wet/dry condition and in-depth absorption of film
            // - film transmissivity used to attenuate radiative emission from solid surface
            //   which assumes same extinction coefficient applies out/in 
            // - operator[](i) is surface temperature of solid region
            scalar qRad = - tabsorptivity[i]*filmRadTransmitted[i] + filmTransmissivity[i]*temissivity[i]*sigma*pow4(operator[](i));

            // Apply additional heat flux at pyrolysis front
            // , if desired
            // Acts to enhance the lateral flame spread
            // , compensating for the 1D heat transfer 
            // , model assumption of pyrolysis region
            scalar extraHeatFlux(0);
            if(qExtraZoneExists && qExtraZone[i] > 0.5)
            {
                extraHeatFlux = (alpha[i]-1.0)*qExtra_;
            }

            // Calculated net heat flux used to set T BC
            nbrTotalFlux[i] = qConv + qRad + extraHeatFlux;

            // Solid temperature update purely based on
            // imposing gradient to satisfy calculated net heat flux
            // applied at solid boundary
            this->refValue()[i] = operator[](i);   
            this->refGrad()[i] = -nbrTotalFlux[i]/K[i];
            this->valueFraction()[i] = 0.0;

            // AK:
            // This code block is purely for setting values read by CUP pyrolysis model
            // complexFuelPyrolysis to be consistent with CUP
            // model as it is currently written while removing the BC complexFuelSuppression
            // which is identical to this BC apart from the following fields being made available
            // to complexFuelPyrolysis.

            qRadiativeTransfer_[i] = filmRadTransmitted[i]; //qRad;
            qConvectiveTransfer_[i] = -qConv;
            qExtraTransfer_[i] = extraHeatFlux;
            wetnessFactorTransfer_[i] = wetnessFraction[i];
            filmAlphaTransfer_[i] = alpha[i];

        }
    }
    else if (mesh.name() == "region0") // gas region 
    {
        
        Tfilm =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "Tsf",
                nbrPatchI,
                true
            );

        alpha =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "alpha",
                nbrPatchI,
                true
            );

        mpp.distribute(Tfilm);
        mpp.distribute(alpha);

        // Twall is blended function of film and solid T
        scalarList Twall(patch().size(), 0.0);

        forAll(*this, i)
        {
            scalar Twet = Tfilm[i]; 
            scalar Tdry = nbrField[i]; // Surface temperature of solid region

            Twall[i] = alpha[i]*(Twet - Tdry) + Tdry;
        }

        // Consistency in T enforced at coupled BC
        this->refValue() = Twall;
        this->refGrad() = 0.0;   
        this->valueFraction() = 1.0;
    }
    else
    {
        FatalErrorIn
        (
            "const filmPyrolysisTemperatureCoupledFFFvPatchScalarField::"
            "updateCoeffs() "
        )
        << "Ill-defined regionName in BC specification, got:  " << mesh.name() 
        << abort(FatalError);
    }

    mixedFvPatchScalarField::updateCoeffs();

}


// These getter functions only used in CUP model 
const scalarField&
filmPyrolysisTemperatureCoupledFFFvPatchScalarField::getRadiativeHeatFlux() const
{
    return qRadiativeTransfer_;
}

const scalarField&
filmPyrolysisTemperatureCoupledFFFvPatchScalarField::getConvectiveHeatFlux() const
{
    return qConvectiveTransfer_;
}

const scalarField&
filmPyrolysisTemperatureCoupledFFFvPatchScalarField::getExtraHeatFlux() const
{
    return qExtraTransfer_;
}

const scalarField&
filmPyrolysisTemperatureCoupledFFFvPatchScalarField::getWetnessFactor() const
{
    return wetnessFactorTransfer_;
}

const scalarField&
filmPyrolysisTemperatureCoupledFFFvPatchScalarField::getFilmAlpha() const
{
    return filmAlphaTransfer_;
}

void filmPyrolysisTemperatureCoupledFFFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>
    (
        "filmRegion",
        "surfaceFilmProperties",
        filmRegionName_
    );
    os.writeEntryIfDifferent<word>
    (
        "pyrolysisRegion",
        "pyrolysisProperties",
        pyrolysisRegionName_
    );
    os.writeKeyword("Tnbr")<< TnbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldRadiativeName")<<
        neighbourFieldRadiativeName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourFieldConvectiveName")<< //Ning
        neighbourFieldConvectiveName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fieldConvectiveName")<< //Ning
        fieldConvectiveName_ << token::END_STATEMENT << nl;
    os.writeKeyword("fieldRadiativeName")<< //Ning
        fieldRadiativeName_ << token::END_STATEMENT << nl;
    os.writeKeyword("K")<< //Ning
        KName_ << token::END_STATEMENT << nl;
    os.writeKeyword("qExtra")<< //Ning
        qExtra_ << token::END_STATEMENT << nl;
    os.writeKeyword("filmDeltaDry")<< // Only used for CUP model
        filmDeltaDry_ << token::END_STATEMENT << nl;
    os.writeKeyword("filmDeltaWet")<< // Only used for CUP model
        filmDeltaWet_ << token::END_STATEMENT << nl;
    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    filmPyrolysisTemperatureCoupledFFFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
