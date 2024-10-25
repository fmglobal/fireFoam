/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::standardRadiation

Description
    Standard radiation model

    Allows for reflection, transmission, and absorption

    Inhereit from this class if special behaviour needed

    In-depth effects follow work presented in:

        "Experimental study on radiation attenuation by a water film"
         D. Brissinger, G. Parent, P. Boulet, JQSRT 145 (2014) 160-168


SourceFiles
    standardRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef standardRadiation_H
#define standardRadiation_H

#include "filmRadiationModel.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                     Class standardRadiation Declaration
\*---------------------------------------------------------------------------*/

class standardRadiation
:
    public filmRadiationModel
{

protected:

        //- Radiative incident flux mapped from  the primary region / [kg/s3]
        volScalarField qinPrimary_;

        //- Diagnostic
        volScalarField blocked_;


private:

    // Private data
        // Model coefficients

            //- Beta coefficient
            scalar beta_;

            //- Bar(kappa) coefficient
            scalar kappaBar_;

        //- Emissivity [-]
        scalar emissivity_;

        //- Minimum fraction of incident radiation NOT blocked[-]
        scalar minNotBlocked_;

        //- MLR (kg/m2/s) at which 60% blockage of radiation expected
        //- larger value leads to weaker blocking as function of MLR
        scalar blockm60_; 

    // Private member functions

        //- No copy construct
        standardRadiation(const standardRadiation&) = delete;

        //- No copy assignment
        void operator=(const standardRadiation&) = delete;
public:

    //- Runtime type information
    TypeName("standardRadiation");


    // Constructors

        //- Construct from surface film model and dictionary
        standardRadiation
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Calling from derived class, pass through the typeName
        standardRadiation
        (
            const word& modelType, 
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~standardRadiation() = default;


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct();

            //- Return the radiation sensible enthalpy source
            virtual tmp<volScalarField> Shs() const;

            //- Return the incident radiation on the film
            virtual tmp<volScalarField> incident() const;

            //- Return the fraction of incident radiation that is not  blocked by fuel vapour above film
            virtual tmp<volScalarField> notBlocked() const; 

            //- Return the radiation emitted from the film
            virtual tmp<volScalarField> emitted() const;

            //- Return the energy absorbed by the film
            virtual tmp<volScalarField> absorbed() const;

            //- Return the energy reflected from the film
            virtual tmp<volScalarField> reflected() const;

            //- Return the energy transmitted through the film to the solid 
            virtual tmp<volScalarField> transmitted() const;

            //- Return the fraction of non-reflected incident energy transmitted through film
            virtual tmp<volScalarField> transmissivity() const;

            //- Return the fraction of energy reflected from the film surface
            virtual scalar albedo() const;

            //- Return the surface emissivity of the film 
            virtual scalar emissivity() const;

            //- Return the extinction coefficient of the liquid material 
            virtual scalar kappaBar() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
