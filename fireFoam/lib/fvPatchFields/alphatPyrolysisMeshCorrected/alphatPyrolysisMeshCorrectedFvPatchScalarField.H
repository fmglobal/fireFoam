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

Class
    alphatPyrolysisMeshCorrectedFvPatchScalarField

Description
    Based on "A convective heat transfer model for LES fire modeling" by
    Ren and Wang PROCI 38(3) 2021

    Modifies SGS thermal diffusivity (alphat) to account for coarse-mesh and blowing
    effects. Modified alphat is used to evaluate convective heat transfer in
    conjunction with resolved thermal boundary layer.

SourceFiles
    alphatPyrolysisMeshCorrectedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphatPyrolysisMeshCorrectedFvPatchScalarField_H
#define alphatPyrolysisMeshCorrectedFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
     Class alphatPyrolysisMeshCorrectedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatPyrolysisMeshCorrectedFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        word   phiName_;
        scalar fuelConversionRatio_;
        scalar m60_;
        scalar Prt_;
        scalar C1_;
        scalar C2_;
        scalar LR_;
        scalar LB_;

public:

    //- Runtime type information
    TypeName("alphatPyrolysisMeshCorrected");


    // Constructors

        //- Construct from patch and internal field
        alphatPyrolysisMeshCorrectedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatPyrolysisMeshCorrectedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given an
        //  alphatPyrolysisMeshCorrectedFvPatchScalarField
        //  onto a new patch
        alphatPyrolysisMeshCorrectedFvPatchScalarField
        (
            const alphatPyrolysisMeshCorrectedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphatPyrolysisMeshCorrectedFvPatchScalarField
        (
            const alphatPyrolysisMeshCorrectedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new alphatPyrolysisMeshCorrectedFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        alphatPyrolysisMeshCorrectedFvPatchScalarField
        (
            const alphatPyrolysisMeshCorrectedFvPatchScalarField&,
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
                new alphatPyrolysisMeshCorrectedFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions

        // Evaluation functions

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::Pstream::commsTypes::blocking
            );

        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
