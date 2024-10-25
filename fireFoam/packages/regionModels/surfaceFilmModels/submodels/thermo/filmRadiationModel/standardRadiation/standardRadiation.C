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

\*---------------------------------------------------------------------------*/

#include "standardRadiation.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardRadiation, 0);

addToRunTimeSelectionTable
(
    filmRadiationModel,
    standardRadiation,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardRadiation::standardRadiation
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    filmRadiationModel(typeName, film, dict),
    qinPrimary_
    (
        IOobject
        (
            "qin", // same name as qin on primary region to enable mapping
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar(dimMass/pow3(dimTime), Zero),
        film.mappedPushedFieldPatchTypes<scalar>()
    ),
    blocked_
    (
        IOobject
        (
            "blocked", 
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    beta_(coeffDict_.lookupOrDefault<scalar>("beta",0.9)),
    kappaBar_(coeffDict_.lookupOrDefault<scalar>("kappaBar",1.0e4)),
    emissivity_(coeffDict_.lookupOrDefault<scalar>("emissivity",0.95)),
    minNotBlocked_(coeffDict_.getCheckOrDefault<scalar>("minNotBlocked",1,scalarMinMax::zero_one())),
    blockm60_(coeffDict_.getOrDefault<scalar>("blockm60",0.01))
    {}


standardRadiation::standardRadiation
(
    const word& modelType, 
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    filmRadiationModel(modelType, film, dict),
    qinPrimary_
    (
        IOobject
        (
            "qin", // same name as qin on primary region to enable mapping
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar(dimMass/pow3(dimTime), Zero),
        film.mappedPushedFieldPatchTypes<scalar>()
    ),
    blocked_
    (
        IOobject
        (
            "blocked", 
            film.time().timeName(),
            film.regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        film.regionMesh(),
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    beta_(coeffDict_.getOrDefault<scalar>("beta",0.9)),
    kappaBar_(coeffDict_.getOrDefault<scalar>("kappaBar",1.0e4)),
    emissivity_(coeffDict_.getOrDefault<scalar>("emissivity",0.95)),
    minNotBlocked_(coeffDict_.getCheckOrDefault<scalar>("minNotBlocked",1,scalarMinMax::zero_one())),
    blockm60_(coeffDict_.getOrDefault<scalar>("blockm60",0.01))
    {}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void standardRadiation::correct()
{
    // Transfer qin from primary region
    qinPrimary_.correctBoundaryConditions();

    blocked_ = 1 - notBlocked();
    blocked_.correctBoundaryConditions();
}


tmp<volScalarField> standardRadiation::Shs() const
{
/*
     Net energy source term in film region
*/
    tmp<volScalarField> tShs
    (
        new volScalarField
        (
            IOobject
            (
                typeName + ":Shs",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimMass/pow3(dimTime), Zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    scalarField& Shs = tShs.ref();

    Shs = absorbed() - emitted();

    return tShs;
}

tmp<volScalarField> standardRadiation::notBlocked() const
{
    tmp<volScalarField> tnotBlocked
    (
        new volScalarField
        (
            IOobject
            (
                typeName + "notBlocked",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimless, Zero)
        )
    );

    scalarField& notBlocked = tnotBlocked.ref();

    scalarField phiFilm =  film().regionMesh().lookupObject<volScalarField>("filmGasMassFlux").primitiveField();

    scalarField mRatio = phiFilm/blockm60_;
    mRatio.clamp_range(ROOTSMALL,10.0);
    scalarField blockFactor = mRatio/(exp(mRatio)-1.0);


    notBlocked = blockFactor + (1-blockFactor)*minNotBlocked_;

    //Info << "phi: " << '\t' << phiFilm << endl;
    //Info << "mRatio: " << '\t' << mRatio << endl;
    //Info << "notBlocked: " << '\t' << notBlocked << endl;
    
    return tnotBlocked;

}


tmp<volScalarField> standardRadiation::incident() const
{
 /*
    Radiation incident on to the film region
 */ 
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
 
    inc =  notBlocked()*qinPrimary_;

    return tinc;
}


tmp<volScalarField> standardRadiation::emitted() const
{
 /*
   Radiation emitted from the film region
 */
    tmp<volScalarField> tem
    (
        new volScalarField
        (
            IOobject
            (
                typeName + "emitted",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimMass/pow3(dimTime), Zero)
        )
    );

	scalarField& em = tem.ref();

    const scalarField& alpha = filmModel_.alpha();

    const scalarField& Tem = qinPrimary_.db().lookupObject<volScalarField>("Tsf"); 

    const scalar sigma = constant::physicoChemical::sigma.value();

    em = alpha*emissivity()*sigma*pow4(Tem);

    return tem;
}


tmp<volScalarField> standardRadiation::absorbed() const
{
/*
    Radiation absorbed by film, not including emited
*/

    tmp<volScalarField> tabs
    (
        new volScalarField
        (
            IOobject
            (
                typeName + "absorbed",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimMass/pow3(dimTime), Zero)
        )
    );

	scalarField& abs = tabs.ref();

    abs = incident() - reflected() - transmitted();

    return tabs;
}


tmp<volScalarField> standardRadiation::reflected() const
{
/*
    Radiation reflected from film
    Alpha accounts for wet/dry 
*/
    tmp<volScalarField> trefl
    (
        new volScalarField
        (
            IOobject
            (
                typeName + "reflected",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimMass/pow3(dimTime), Zero)
        )
    );

	scalarField& refl = trefl.ref();

    const scalarField& alpha = filmModel_.alpha();

    refl = alpha*incident()*albedo(); 

    return trefl;
}


// Energy transmitted through film to solid
tmp<volScalarField> standardRadiation::transmitted() const
{
/*
    Radiation transmitted through film to the solid region
    Reflected component and transmissivity calculations 
    account for wet/dry condition    
*/
    tmp<volScalarField> ttran
    (
        new volScalarField
        (
            IOobject
            (
                typeName + "transmitted",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimMass/pow3(dimTime), Zero)
        )
    );
	
    scalarField& tran = ttran.ref();

    tran = ( incident() - reflected() ) * transmissivity();

    return ttran;
}

// Fraction of non-reflected incident energy transmitted to solid
tmp<volScalarField> standardRadiation::transmissivity() const
{
/*
    Calculate transmissivity of film
    Accounts for extinction coefficient and film depth
    Applies alpha to account for wet/dry threshold
*/ 

    tmp<volScalarField> ttau
    (
        new volScalarField
        (
            IOobject
            (
                typeName + "transmissivity",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

	const scalarField& delta = filmModel_.delta();

    const scalarField& alpha = filmModel_.alpha();

	scalarField& tau = ttau.ref();

    tau = min(max(exp(-kappaBar()*delta),SMALL),1.0);

    /*
     alpha = 0 --> dry, tau = 1, all energy transmitted
     alpha = 1 --> wet, tau = tau, depth-dependent trans. calc unchanged
     0<alpha<1 --> partial wet, blended dry and wet calc.
    */
    tau = alpha * ( tau - 1 ) + 1;

    return ttau;
}


// Fraction of incident energy reflected from surface of film 
scalar standardRadiation::albedo() const
{
    return 1 - beta_;
}

// Surface emissivity of the film
scalar standardRadiation::emissivity() const
{
    return emissivity_;
}

// Extinction coefficient units m-1
scalar standardRadiation::kappaBar() const
{
    return kappaBar_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
