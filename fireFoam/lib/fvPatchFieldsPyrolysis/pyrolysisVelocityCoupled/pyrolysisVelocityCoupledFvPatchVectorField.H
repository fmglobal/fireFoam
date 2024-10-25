/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2010 OpenCFD Ltd.
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

Class
    Foam::pyrolysisVelocityCoupledFvPatchVectorField

Description
    Based on mappedFlowRateFvPatchVectorField with modifications for HOC
    correction.

    Coupled the mass loss from a pyrolysis region to the velocity inlet into
    the primary (gas) region.

    The mass flux phi is converted into a velocity through the relationship
    U = (mdot/rho/dA)

    The velocity must be corrected based on the knowledge of the chemical 
    difference between the notional gasification products (pyrolysate) and
    the surrogate gas phase species used to represent the pyrolysate.

    For example, typical cellulosic fuels have a heat of combustion (HOC) of
    15-20 MJ/kg, while typical combustion fuels used to represent these 
    complex mixtures have HOC values of 40-50 MJ/kg.

    To predict the correct fire size (energy content of the gasses introduced
    into the gas-phase domain), the velocity BC must be corrected:

    U_BC = U ( HOC_pyr / HOC_g ), 

    Where HOC_pyr is the HOC of the pyrolosate mixture (user input) and 
    HOC_g is the HOC of the gas species, determined by looking-up from
    the thermo object that handles the gas-phase thermo-chemistry.

    Two further corrections are supported:

    1) charring_ (required input): corrects the HOC_pyr value for correcting 
    difference between measured pyrolysate HOC (includes char oxidation)
    and the volatile component only, excluding char oxidization which is more
    relavent for time-scales of typical fire simulation (persistent char layer
    is slow to oxidize).

    2) complexFuel_ (default false): corrects HOC_pyr for complex commodities that
    have time-varying pyrolysate composition. 

SourceFiles
    pyrolysisVelocityCoupledFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pyrolysisVelocityCoupledFvPatchVectorField_H
#define pyrolysisVelocityCoupledFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class flowRateInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class pyrolysisVelocityCoupledFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Name of the neighbor flux setting the inlet mass flux
        word nbrPhiName_;

        //- Name of the local mass flux
        word phiName_;

        //- Name of the density field used to normalize the mass flux
        word rhoName_;

        //- Name of the local mass flux
        word virginName_;

        //- Name of the density field used to normalize the mass flux
        word charName_;

        //- Heat of combustion of virgin pyrolysate
        const scalar hocVirgin_;

        //- Head of combustion of char pyrolysate
        const scalar hocChar_;

        //- Include correction for char oxidation
        const bool charring_;

        //- Flag to identify CUP/Class-3 regions
        const bool complexFuel_;


public:

   //- Runtime type information
   TypeName("pyrolysisVelocityCoupled");


   // Constructors

        //- Construct from patch and internal field
        pyrolysisVelocityCoupledFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pyrolysisVelocityCoupledFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  pyrolysisVelocityCoupledFvPatchVectorField
        //  onto a new patch
        pyrolysisVelocityCoupledFvPatchVectorField
        (
            const pyrolysisVelocityCoupledFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pyrolysisVelocityCoupledFvPatchVectorField
        (
            const pyrolysisVelocityCoupledFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new pyrolysisVelocityCoupledFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pyrolysisVelocityCoupledFvPatchVectorField
        (
            const pyrolysisVelocityCoupledFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new pyrolysisVelocityCoupledFvPatchVectorField(*this, iF)
            );
        }


    // Member functions


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
