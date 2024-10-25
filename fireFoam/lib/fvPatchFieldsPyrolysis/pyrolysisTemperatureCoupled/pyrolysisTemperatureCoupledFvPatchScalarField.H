/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::
    compressible::
    pyrolysisTemperatureCoupledFvPatchScalarField

Description
    Temperture BC for coupling pyrolysis region with primary region.

    Consider the gas-solid interface, the assumption is that the 
    net heat transfer from the gas-phase into the interface is
    matched by the heat transfer from the interface into the solid.

    i.e., net heat flux into solid (convective + radiative) equals
    heat flux conducted into solid.

    This is achevied by setting the T BC to result in a conduction
    term equal to the calculated net flux on the solid side, on the
    gas side the temperature is matched to the same value.

SourceFiles
    pyrolysisTemperatureCoupledFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef pyrolysisTemperatureCoupledFvPatchScalarField_H
#define pyrolysisTemperatureCoupledFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "pyrolysisModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
   Class pyrolysisTemperatureCoupledFvPatchScalarField
   Declaration
\*---------------------------------------------------------------------------*/

class pyrolysisTemperatureCoupledFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
public:

    typedef Foam::regionModels::pyrolysisModels::pyrolysisModel
        pyrolysisModelType;

private:

    // Private data

       //- Name of pyrolysis region
        const word pyrolysisRegionName_;

        //- Name of field on the neighbour region
        const word TnbrName_;

         //- Name of the radiative heat flux on the neighbout region
        const word neighbourFieldRadiativeName_;
        const word neighbourFieldConvectiveName_;

        //- Name of thermal conductivity field
        const word KName_;

public:

    //- Runtime type information
    TypeName("compressible::pyrolysisTemperatureCoupled");


    // Constructors

        //- Construct from patch and internal field
        pyrolysisTemperatureCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pyrolysisTemperatureCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        pyrolysisTemperatureCoupledFvPatchScalarField
        (
            const
            pyrolysisTemperatureCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new pyrolysisTemperatureCoupledFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        pyrolysisTemperatureCoupledFvPatchScalarField
        (
            const pyrolysisTemperatureCoupledFvPatchScalarField&,
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
                new pyrolysisTemperatureCoupledFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        virtual void autoMap(const fvPatchFieldMapper&);

        virtual void rmap(const fvPatchScalarField& psf, const labelList& addr);

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
