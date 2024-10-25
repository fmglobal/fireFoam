/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::radiation::RVFGasSootEmission

Description
    RVFGasSootEmission radiation model coupled with the RVF flame extinction model

    The local flame radiant fraction is calculated accounting for gas and soot
    emission. The effective soot volume fraction is simply modeled using the
    adiabatic flame temperature. Flame thickness is calculated in the RVF model.

SourceFiles
    RVFGasSootEmission.C

\*---------------------------------------------------------------------------*/

#ifndef RVFGasSootEmission_H
#define RVFGasSootEmission_H

#include "interpolationLookUpTable.H"
#include "absorptionEmissionModel.H"
#include "HashTable.H"
#include "absorptionCoeffs.H"
#include "basicThermo.H"

//#include "OpenFOAM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                 Class RVFGasSootEmission Declaration
\*---------------------------------------------------------------------------*/

class RVFGasSootEmission
:
    public absorptionEmissionModel
{
public:

    // Public data

        // Maximum number of species considered for absorptivity
        static const int nSpecies_ = 5;

        // Absorption Coefficients
        absorptionCoeffs coeffs_[nSpecies_];


private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Emission constant coefficient
        const scalar EhrrCoeff_;

        //- Switch for scaling radiant fraction
        Switch radScaling;

        //- flag to indicate whether patch2 is a CUP region
        const Switch patch2isCUP_;

        //- Radiant fraction for fuel from patch 1 and 2
        scalar Ehrr1_;
        scalar Ehrr2_;

        //- patch names for patches 1 and 2
        wordList patchName1_;
        wordList patchName2_;

        //- CUP region: radiant fractions for corrugated Carbon and polystyrene
        const scalar Ehrr2CC_;
        const scalar Ehrr2PS_;

        //- min and max specified fuel radiant fractions
        const scalar EhrrMin_;
        const scalar EhrrMax_;

	//- for CUP: time to start using local radiant fractions
	const scalar tstartXrLocal_;

public:

    //- Runtime type information
    TypeName("RVFGasSootEmission");


    // Constructors

        //- Construct from components
        RVFGasSootEmission
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~RVFGasSootEmission();


    // Member Operators

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> aCont(const label bandI = 0) const;


            // Emission coefficient

                //- Emission coefficient for continuous phase
                tmp<volScalarField> eCont(const label bandI = 0) const;


            // Emission contribution

                //- Emission contribution for continuous phase
                tmp<volScalarField> ECont(const label bandI = 0) const;


    // Member Functions

        inline bool isGrey() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
