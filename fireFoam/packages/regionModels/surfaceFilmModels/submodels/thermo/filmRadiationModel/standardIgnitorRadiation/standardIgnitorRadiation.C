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

\*---------------------------------------------------------------------------*/

#define DEBUG(x) {                                              \
        std::streamsize p = std::cout.precision();              \
        std::ios::fmtflags myFlags;                             \
        myFlags = cout.flags();                                 \
        std::cout.precision(10);                                \
        std::cout.setf(std::ios::fixed,std::ios::floatfield);   \
        std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "; \
        std::cout << "p" << Pstream::myProcNo();                \
        std::cout << " " << #x " = " << x << std::endl;         \
        std::cout.precision(p);                                 \
        std::cout.flags(myFlags);                               \
    }
#define TRACE(s) {                                              \
        std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "; \
        std::cout << "p" << Pstream::myProcNo();                \
        std::cout << " " << #s << std::endl;                    \
        s;                                                      \
    }

#include "standardIgnitorRadiation.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardIgnitorRadiation, 0);

addToRunTimeSelectionTable
(
    filmRadiationModel,
    standardIgnitorRadiation,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardIgnitorRadiation::standardIgnitorRadiation
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    standardRadiation(typeName, film, dict),
    qrConst_
    (
        IOobject
        (
            typeName + ":qrConst",
            film.time().timeName(),
            film.regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh()
    ),
    mask_
    (
        IOobject
        (
            typeName + ":mask",
            film.time().timeName(),
            film.regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar("one", dimless, 1.0)
    ),
    timeStart_(readScalar(coeffDict_.lookup("timeStart"))),
    duration_(readScalar(coeffDict_.lookup("duration")))
{
    mask_ = pos0(mask_ - 0.5); // kvm
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardIgnitorRadiation::~standardIgnitorRadiation()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void standardIgnitorRadiation::correct() 
{
   standardRadiation::correct(); 
}

// Incident radiation determined from user-specified BCs
tmp<volScalarField> standardIgnitorRadiation::incident() const
{
    tmp<volScalarField> tinc
    (
        new volScalarField
        (
            IOobject
            (
                typeName + "incident",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimMass/pow3(dimTime), Zero)
        )
    );
	scalarField& inc = tinc.ref();

    const scalar time = film().time().value();

    if ( active(time) )
    {
        // Mask:
        //       0 -> do not apply fixed rad from BC, use qinP from RTE solve
        //       1 -> use fixed rad from BC, overide qinP from RTE solve
        inc = notBlocked()*(mask_*qrConst_ + (1-mask_)*qinPrimary_);
    }
    else
    {
        inc = notBlocked()*qinPrimary_;
    }


    return tinc;
}

bool standardIgnitorRadiation::active(scalar time) const
{

    if ( (time >= timeStart_) && (time <= timeStart_ + duration_ ) )
    {
        return true;
    }
    else
    {
        return false;
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
