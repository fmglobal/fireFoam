/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::standardIgnitorRadiation

Description
    Combination of standardRadiation and standardIgnitorRadiation behaviour

    If time is between user defiend start and end points:

        - where mask_ = 0, set incident() to be qinPrimary_
        - where mask_ = 1, set incident() to be qrConst_

    Else:

        - set incident() to be qinPrimary_

    This model allows for user to give an ignition 'hot-spot' to initiate
    the interesting physics, but still take the qin value form the RTE
    solution at times and locations other than those specified by the user.

    Note - the user defined hotspot over-writes the qin value from the RTE,
    it does not add to it (prevent unrealistically large incident() values).

SourceFiles
    standardIgnitorRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef standardIgnitorRadiation_H
#define standardIgnitorRadiation_H

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
                      Class standardIgnitorRadiation Declaration
\*---------------------------------------------------------------------------*/

class standardIgnitorRadiation
:
    public standardRadiation 
{
private:

    // Private data

        //- Constant radiative flux [kg/s3]
        volScalarField qrConst_;

        //- Radiation mask
        volScalarField mask_;

        //- Time start [s]
        const scalar timeStart_;

        //- Duration [s]
        const scalar duration_;

    // Private member functions

        //- No copy construct
        standardIgnitorRadiation(const standardIgnitorRadiation&) = delete;

        //- No copy assignment
        void operator=(const standardIgnitorRadiation&) = delete;


public:

    //- Runtime type information
    TypeName("standardIgnitorRadiation");

    // Constructors

        //- Construct from surface film model and dictionary
        standardIgnitorRadiation
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~standardIgnitorRadiation();

    // Member Functions

    virtual void correct();

    //- Return the incident radiation on the film
    virtual tmp<volScalarField> incident() const;

    bool active(scalar) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
