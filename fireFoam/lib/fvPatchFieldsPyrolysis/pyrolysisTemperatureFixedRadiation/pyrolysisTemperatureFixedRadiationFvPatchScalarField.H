/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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

Class
    Foam::pyrolysisTemperatureFixedRadiationFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    Boundary condition for imposing radiative heat flux for solid regions.
    Used to emulate a fixed incident radiative heat flux on a wall.

    the gradient heat flux is calculated as :

    qr = ( absorptivity*qrIncident - emissivity*sigma_*T^4)

    where:

    emissivity is the emissivity of the solid.
    absorptivity is the absorptivity of the solid.
    qrIncident is the specified fixed incident radiation.

    Example usage:

    wall
    {
        type            pyrolysisTemperatureFixedRadiation;
        qrIncident      uniform 500;
        kappa           solidThermo;
        KappaName       none;
    }

    kappa:
    - 'lookup' :    lookup volScalarField (or volSymmTensorField) with name
    - 'solidThermo' : use solidThermo kappa()

    emissivity:
    - 'lookup' :     lookup volScalarField emissivity
    - 'localSolidRadiation': Look up for local solidRadiation


SourceFiles
    pyrolysisTemperatureFixedRadiationFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef pyrolysisTemperatureFixedRadiationFvPatchScalarField_H
#define pyrolysisTemperatureFixedRadiationFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"
#include "temperatureCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class pyrolysisTemperatureFixedRadiationFvPatchScalarField declaration
\*---------------------------------------------------------------------------*/

class pyrolysisTemperatureFixedRadiationFvPatchScalarField
:
    public fixedGradientFvPatchScalarField,
    public temperatureCoupledBase
{
    // Private data

        //- Incident radiative heat flux
        scalarField qrIncident_;


public:

    //- Runtime type information
    TypeName("pyrolysisTemperatureFixedRadiation");


    // Constructors

        //- Construct from patch and internal field
        pyrolysisTemperatureFixedRadiationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pyrolysisTemperatureFixedRadiationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        pyrolysisTemperatureFixedRadiationFvPatchScalarField
        (
            const
            pyrolysisTemperatureFixedRadiationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pyrolysisTemperatureFixedRadiationFvPatchScalarField
        (
            const pyrolysisTemperatureFixedRadiationFvPatchScalarField&
        );


        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new pyrolysisTemperatureFixedRadiationFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        pyrolysisTemperatureFixedRadiationFvPatchScalarField
        (
            const pyrolysisTemperatureFixedRadiationFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new pyrolysisTemperatureFixedRadiationFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions


     // Mapping functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(const fvPatchFieldMapper&);

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchScalarField&,
            const labelList&
        );

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
