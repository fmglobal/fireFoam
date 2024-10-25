/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015 OpenCFD Ltd.
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

#include "greyMeanSolidAbsorptionEmissionFF.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(greyMeanSolidAbsorptionEmissionFF, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            greyMeanSolidAbsorptionEmissionFF,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::radiation::
greyMeanSolidAbsorptionEmissionFF::X(const word specie) const
{
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();

    tmp<scalarField> tXj(new scalarField(T.primitiveField().size(), Zero));
    scalarField& Xj = tXj.ref();

    tmp<scalarField> tRhoInv(new scalarField(T.primitiveField().size(), Zero));
    scalarField& rhoInv = tRhoInv.ref();

    forAll(mixture_.Y(), specieI)
    {
        const scalarField& Yi = mixture_.Y()[specieI];

        forAll(rhoInv, iCell)
        {
            rhoInv[iCell] +=
                Yi[iCell]/mixture_.rho(specieI, p[iCell], T[iCell]);
        }
    }
    const scalarField& Yj = mixture_.Y(specie);
    const label mySpecieI = mixture_.species().find(specie);
    forAll(Xj, iCell)
    {
        Xj[iCell] = Yj[iCell]/mixture_.rho(mySpecieI, p[iCell], T[iCell]);
    }

    return (Xj/rhoInv);
}

// alex
void Foam::radiation::
greyMeanSolidAbsorptionEmissionFF::propertyPolynomial
(
    const word specie,
    const label propertyId
) const
{
    const volScalarField& T = thermo_.T();

    //volScalarField& temporaryPropertyRef = temporaryProperty_;
    const label specieI = mixture_.species()[specie];
    forAll(temporaryProperty_, iCell)
    {
        temporaryProperty_[iCell] = evaluatePolynomial(specieI, propertyId, T[iCell]);
    }

}

// alex
Foam::scalar Foam::radiation::
greyMeanSolidAbsorptionEmissionFF::evaluatePolynomial
(
    const label specie,
    const label property,
    const scalar Ti
) const
{
    const FixedList<scalar,7>& polyCoeffs = polyData_[specie][property];

    const scalar T(max(polyMinTemperature_,min(polyMaxTemperature_,Ti)));

    return polyCoeffs[0]+polyCoeffs[1]*T+polyCoeffs[2]*pow(T,2)
          +polyCoeffs[3]*pow(T,3)+polyCoeffs[4]*pow(T,4)
          +polyCoeffs[5]*pow(T,5)+polyCoeffs[6]*pow(T,6);
} 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyMeanSolidAbsorptionEmissionFF::
greyMeanSolidAbsorptionEmissionFF
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.optionalSubDict(typeName + "Coeffs"))),
    thermo_(mesh.lookupObject<solidThermo>(basicThermo::dictName)),
    speciesNames_(0),
    mixture_(dynamic_cast<const basicSpecieMixture&>(thermo_)),
 // alex
    constData_(mixture_.Y().size()),
    usePolynomial_(mixture_.Y().size()),
    polyData_(mixture_.Y().size()),
    polyMaxTemperature_(coeffsDict_.getOrDefault<scalar>("polyMaxTemperature",3000)),
    polyMinTemperature_(coeffsDict_.getOrDefault<scalar>("polyMaxTemperature",300)),
    temporaryProperty_
    (
        IOobject
        (
            "temporaryRadiationProperty",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero),
        extrapolatedCalculatedFvPatchVectorField::typeName
    )
{
    if (!isA<basicSpecieMixture>(thermo_))
    {
        FatalErrorInFunction
            << "Model requires a multi-component thermo package"
            << abort(FatalError);
    }

    label specie = 0;
    const dictionary& functionDicts = dict.optionalSubDict(typeName + "Coeffs");

    Info << nl <<  "Setting radiative properties for region " << mesh.name() << nl << endl;

    // Iterate over each species in Coeffs dict
    for (const entry& eSpecie : functionDicts)
    {
        if (!eSpecie.isDict())  // safety
        {
            continue;
        }

        const word& specieKey = eSpecie.keyword();
        const dictionary& dSpecie = eSpecie.dict();

        if (!mixture_.contains(specieKey))
        {
            WarningInFunction
                << " specie: " << specieKey << " is not found in the solid mixture"
                << nl
                << " specie is the mixture are:" << mixture_.species() << nl
                << nl << endl;
        }

        Info << "Species: " << specieKey << endl;

        speciesNames_.insert(specieKey, specie);

        // Iterate over each property in specie dict
        for (const entry& eProperty : dSpecie)
        {

            // Set enum value for property
            const word& propertyKey = eProperty.keyword();
            
            Info << tab << "Property: " << propertyKey << endl;

            label property(Zero);
            const word& abs("absorptivity");
            const word& em("emissivity");

            if(propertyKey == abs)
            {
                property = absorptivity;
            }
            else if (propertyKey == em)
            {
                property = emissivity;
            }
            else
            {
                FatalErrorInFunction
                    << " unknown property " << propertyKey << " specified in " << nl
                    << " dictionary " << functionDicts.name() << abort(FatalError);
            }

            // Don't use polynomial evaluation unless needed
            usePolynomial_[specie][property] = false;
            
            // Support old method where scalar value assigned 
            if (!eProperty.isDict())
            {
                dSpecie.readEntry<scalar>(propertyKey,constData_[specie][property]);

                Info << tab << tab << "Method: " << "DEPRECATED" << nl
                     << tab << tab << "Value: " << constData_[specie][property] << endl;

                WarningInFunction
                    << " specie: " << specieKey << " does not have a dictionary "
                    << "for property " << propertyKey << nl 
                    << "Reverting to deprecated input format assuming species dict contains "
                    << "fixed scalar value for this property" << nl
                    << "Update to new input method as follows..." << nl
                    << "absorptivity" << nl
                    << "{" << nl
                    << "    method  <fixedValue, fixedTemperature, polynomial>;" << nl
                    << "    fixedValue  <value>;  // Needed if fixedValue selected" << nl
                    << "    fixedTemperature <Temperature>; // If fixedTemperature selected" << nl
                    << "    polynomialCoeffs <[a b c d e f g]; // Needed if fixedValue is not selected" << nl
                    << "}"
                    << endl;

                continue;
            }

            const dictionary& dProperty = eProperty.dict();

            word methodKey;
            label method(fixedValue);
            FixedList<scalar,7> polyCoeffs;

            dProperty.readEntry("method",methodKey);

            if(methodKey=="fixedValue")
            {
                method = fixedValue;
            }
            else if (methodKey=="fixedTemperature")
            {
                method = fixedTemperature;
            }
            else if (methodKey == "polynomial")
            {
                if (property == absorptivity)
                {
                    // Polynomial evaluation uses variable solid surface temperature
                    // Does not make sense for absorptivity, better to presume 
                    // fixed temperature of radiative source (e.g. flame, lamp)
                    WarningInFunction << "Polynomial evaluation for absorptivity is not currently supported. Reverting to fixedTemperature method." << endl;
                    methodKey = "fixedTemperature";
                    method = fixedTemperature;
                }
                else if (property == emissivity)
                {
                    method = polynomial;
                    usePolynomial_[specie][property] = true;
                }
                else
                {
                    FatalErrorInFunction << "Undefined property" << abort(FatalError);
                }
            }
            else
            {
                FatalErrorInFunction
                    << " unknown method " << methodKey << " for property " << propertyKey
                    << " of species " << specieKey << " specified in " <<  " dictionary " << functionDicts.name()
                    << abort(FatalError);     
            }

            // Read polynomial data if it's relevant
            if (method == polynomial || method == fixedTemperature)
            {
                dProperty.readEntry<FixedList<scalar,7>>("polynomialCoeffs",polyCoeffs);

                forAll(polyCoeffs,c)
                {
                    polyData_[specie][property][c] = polyCoeffs[c];
                }    
            }

            switch(method)
            {
                case fixedValue:
                    dProperty.readEntry("fixedValue",constData_[specie][property]);
                    Info << tab << tab << "Method: " << methodKey << nl
                         << tab << tab << "Value: " << constData_[specie][property] << endl;
                    break;
                case polynomial:
                    // Values assigned on the fly based on solid temperature
                    usePolynomial_[specie][property] = true;
                    Info << tab << tab << "Method: " << methodKey << nl
                         << tab << tab << "Value: " << "dynamic" << endl;
                    break;
                case fixedTemperature:
                    scalar propertyEvaluationTemperature;
                    dProperty.readEntry("fixedTemperature",propertyEvaluationTemperature);
                    constData_[specie][property] = evaluatePolynomial(specie,property,propertyEvaluationTemperature);
                    Info << tab << tab << "Method: " << methodKey << nl
                         << tab << tab << "Characteristic temperature: " << propertyEvaluationTemperature << " (K)" << nl
                         << tab << tab << "Value: " << constData_[specie][property] << endl;
                    break;
                default:
                    FatalErrorInFunction
                        << " unknown method attempted for property evaluation" << abort(FatalError);
                    break;
            }
        }

        ++specie;
    }
    Info << nl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanSolidAbsorptionEmissionFF::
calc(const label propertyId) const
{
    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar(dimless/dimLength, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );

    scalarField& a = ta.ref().primitiveFieldRef();

    forAllConstIters(speciesNames_, iter)
    {
        if (mixture_.contains(iter.key()))
        {
            if (usePolynomial_[iter()][propertyId])
            {
                // AK - this method is ugly, but it works
                //    - the cleaner approach commented out
                //      below is buggy, see note
                propertyPolynomial(iter.key(),propertyId); // sets temporaryProperty_
                a += temporaryProperty_*X(iter.key());

                // AK - in this version, propertyPolynomial returned  
                //      a tmp<scalarField>, in the same way as for X(iter.key())
                //    - for an unknown reason, the size of the return value object
                //      was wrong, this lead to junk values being assigned without
                //      warning or error.
                //a += propertyPolynomial(iter.key(),propertyId)*X(iter.key());
            }
            else
            {
                a += constData_[iter()][propertyId]*X(iter.key());
            }
        }
    }

    ta.ref().correctBoundaryConditions();
    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanSolidAbsorptionEmissionFF::eCont
(
    const label bandI
) const
{
   return calc(emissivity);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyMeanSolidAbsorptionEmissionFF::aCont
(
    const label bandI
) const
{
   return calc(absorptivity);
}

// ************************************************************************* //
