/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::rampingRadiation

Description
    Ramping radiation model. Starting at timeStart_ and running for duration_
    applies qrConst_ to the film where exposed by mask_ values

     

    qrConst_ and mask_ are set in the boundary conditions



SourceFiles
    rampingRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef rampingRadiation_H
#define rampingRadiation_H

#include "standardRadiation.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class rampingRadiation Declaration
\*---------------------------------------------------------------------------*/

class rampingRadiation
:
    public standardRadiation 
{
private:

    // Private data

        //- Constant radiative flux [kg/s3]
        mutable volScalarField qrConst_;

        //- Radiation mask
        volScalarField mask_;

        //- The rate at which the incident radiation is increasing [kW/m2/s]  
        scalar rampRate_;

        //- Time start [s]
        const scalar timeStart_;

        //- Duration [s]
        const scalar duration_;


    // Private member functions

        //- Disallow default bitwise copy construct
        rampingRadiation(const rampingRadiation&);

        //- Disallow default bitwise assignment
        void operator=(const rampingRadiation&);


public:

    //- Runtime type information
    TypeName("rampingRadiation");


    // Constructors

        //- Construct from surface film model and dictionary
        rampingRadiation
        (
            surfaceFilmRegionModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~rampingRadiation();


    // Member Functions

    //- Return the radiation sensible enthalpy source
    virtual tmp<volScalarField> incident() const;

    bool rampTime(scalar) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
