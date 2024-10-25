/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "liquidFilmThermo.H"
#include "demandDrivenData.H"
#include "thermoSingleLayer.H"
#include "SLGThermo.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(liquidFilmThermo, 0);

addToRunTimeSelectionTable
(
    filmThermoModel,
    liquidFilmThermo,
    dictionary
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const thermoSingleLayer& liquidFilmThermo::thermoFilm() const
{
    if (!isA<thermoSingleLayer>(filmModel_))
    {
        FatalErrorInFunction
            << "Thermo model requires a " << thermoSingleLayer::typeName
            << " film to supply the pressure and temperature, but "
            << filmModel_.type() << " film model selected.  "
            << "Use the 'useReferenceValues' flag to employ reference "
            << "pressure and temperature" << exit(FatalError);
    }

    return refCast<const thermoSingleLayer>(filmModel_);
}


void liquidFilmThermo::initLiquid(const dictionary& dict)
{
    if (liquidPtr_ != nullptr)
    {
        return;
    }

    dict.readEntry("liquid", name_);

    const SLGThermo* thermoPtr =
        filmModel_.primaryMesh().findObject<SLGThermo>("SLGThermo");

    if (thermoPtr)
    {
        // Retrieve from film thermo
        ownLiquid_ = false;

        const SLGThermo& thermo = *thermoPtr;

        const label id = thermo.liquidId(name_);

        liquidPtr_ = &thermo.liquids().properties()[id];
    }
    else
    {
        // New liquid create
        ownLiquid_ = true;

        liquidPtr_ =
            liquidProperties::New(dict.optionalSubDict(name_ + "Coeffs")).ptr();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

liquidFilmThermo::liquidFilmThermo
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    filmThermoModel(typeName, film, dict),
    name_("unknown_liquid"),
    liquidPtr_(nullptr),
    ownLiquid_(false),
    useReferenceValues_(coeffDict_.get<bool>("useReferenceValues")),
    useSurfaceTemperature_(coeffDict_.getOrDefault<bool>("useSurfaceTemperature",false)),
    muScale_(coeffDict_.getOrDefault<scalar>("muScale",Zero)),
    muCritT_(coeffDict_.getOrDefault<scalar>("muCritT",473.0)),
    pRef_(coeffDict_.getOrDefault<scalar>("pRef",101325.0)),
    Tref_(coeffDict_.getOrDefault<scalar>("Tref",298.15)), 
    Tmin_(coeffDict_.getOrDefault<scalar>("Tmin",-VGREAT)),
    Tmax_(coeffDict_.getOrDefault<scalar>("Tmax",VGREAT)) 
{
    initLiquid(coeffDict_);

    Info<< "    film thermo reference temmperature  (K): " << Tref_ << endl; 
    Info<< "    film thermo reference pressure  (Pa): " << pRef_ << endl; 

    Info<< "    limiting minimum film temperature for property evaluaiton to " << Tmin_ << endl; 
    Info<< "    limiting maximum film temperature for property evaluation to " << Tmax_ << endl;

    if (useSurfaceTemperature_)
    {
        Info<< "\tUsing surface temperature for transport property evaluation" << endl; 
    }
  
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

liquidFilmThermo::~liquidFilmThermo()
{
    if (ownLiquid_)
    {
        deleteDemandDrivenData(liquidPtr_);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const word& liquidFilmThermo::name() const
{
    return name_;
}

// Alex
scalar liquidFilmThermo::Tref() const
{
    return Tref_;
}

// Alex
scalar liquidFilmThermo::pRef() const
{
    return pRef_;
}

scalar liquidFilmThermo::rho
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->rho(p, T);
}

// Alex
scalar liquidFilmThermo::muPrefactor
(
    const scalar T
) const
{

    const scalar Tdiff(max(0,T-muCritT_));

    return exp(-muScale_*Tdiff); 
}

scalar liquidFilmThermo::mu
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->mu(p, T);
}


scalar liquidFilmThermo::sigma
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->sigma(p, T);
}


scalar liquidFilmThermo::Cp
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->Cp(p, T);
}


scalar liquidFilmThermo::kappa
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->kappa(p, T);
}


scalar liquidFilmThermo::D
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->D(p, T);
}


scalar liquidFilmThermo::hl
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->hl(p, T);
}


scalar liquidFilmThermo::pv
(
    const scalar p,
    const scalar T
) const
{
    return liquidPtr_->pv(p, T);
}


scalar liquidFilmThermo::W() const
{
    return liquidPtr_->W();
}


scalar liquidFilmThermo::Tb(const scalar p) const
{
    return liquidPtr_->pvInvert(p);
}

scalar liquidFilmThermo::hs
(
    const scalar p,
    const scalar T,
    const bool withRespectToRef = true
) const
{

    scalar hs(liquidPtr_->h(p,T));

    if (withRespectToRef)
    {
        hs -= liquidPtr_->h(this->pRef(),this->Tref());
    }

    return hs;
}

tmp<volScalarField> liquidFilmThermo::rho() const
{
    tmp<volScalarField> trho
    (
        new volScalarField
        (
            IOobject
            (
                type() + ":rho",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimDensity, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );

    scalarField& rho = trho.ref().primitiveFieldRef();

    if (useReferenceValues_)
    {
        rho = this->rho(pRef_, Tref_);
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& p = film.pPrimary();

        forAll(rho, celli)
        {
            const scalar Ts = max(min(T[celli],Tmax_),Tmin_); // kvm
            rho[celli] = this->rho(p[celli], Ts);
        }
    }

    trho.ref().correctBoundaryConditions();

    return trho;
}


tmp<volScalarField> liquidFilmThermo::mu() const
{
    tmp<volScalarField> tmu
    (
        new volScalarField
        (
            IOobject
            (
                type() + ":mu",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimPressure*dimTime, Zero),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    scalarField& mu = tmu.ref().primitiveFieldRef();

    if (useReferenceValues_)
    {
        mu = this->mu(pRef_, Tref_);
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& Ts = film.Ts();
        const volScalarField& p = film.pPrimary();

        scalar Tv(Zero);

        forAll(mu, celli)
        {
            Tv = useSurfaceTemperature_ ? Ts[celli] : T[celli];
            Tv = max(min(Tv,Tmax_),Tmin_); // kvm
            mu[celli] = muPrefactor(Tv)*this->mu(p[celli], Tv); // kvm
        }
    }

    tmu.ref().correctBoundaryConditions();

    return tmu;
}


tmp<volScalarField> liquidFilmThermo::sigma() const
{
    tmp<volScalarField> tsigma
    (
        new volScalarField
        (
            IOobject
            (
                type() + ":sigma",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimMass/sqr(dimTime), Zero),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    scalarField& sigma = tsigma.ref().primitiveFieldRef();

    if (useReferenceValues_)
    {
        sigma = this->sigma(pRef_, Tref_);
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& Ts = film.Ts();
        const volScalarField& p = film.pPrimary();
    
        scalar Tv(Zero);

        forAll(sigma, celli)
        {
            Tv = useSurfaceTemperature_ ? Ts[celli] : T[celli]; 
            Tv = max(min(Tv,Tmax_),Tmin_); // kvm
            sigma[celli] = this->sigma(p[celli], Tv); // kvm
        }
    }

    tsigma.ref().correctBoundaryConditions();

    return tsigma;
}


tmp<volScalarField> liquidFilmThermo::Cp() const
{
    tmp<volScalarField> tCp
    (
        new volScalarField
        (
            IOobject
            (
                type() + ":Cp",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimEnergy/dimMass/dimTemperature, Zero),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    scalarField& Cp = tCp.ref().primitiveFieldRef();

    if (useReferenceValues_)
    {
        Cp = this->Cp(pRef_, Tref_);
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& p = film.pPrimary();

        forAll(Cp, celli)
        {
            const scalar Ts = max(min(T[celli],Tmax_),Tmin_); // kvm
            Cp[celli] = this->Cp(p[celli], Ts); // kvm
        }
    }

    tCp.ref().correctBoundaryConditions();

    return tCp;
}


tmp<volScalarField> liquidFilmThermo::kappa() const
{
    tmp<volScalarField> tkappa
    (
        new volScalarField
        (
            IOobject
            (
                type() + ":kappa",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimPower/dimLength/dimTemperature, Zero),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    scalarField& kappa = tkappa.ref().primitiveFieldRef();

    if (useReferenceValues_)
    {
        kappa = this->kappa(pRef_, Tref_);
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();

        const volScalarField& T = film.T();
        const volScalarField& p = film.pPrimary();

        forAll(kappa, celli)
        {
            const scalar Ts = max(min(T[celli],Tmax_),Tmin_); // kvm
            kappa[celli] = this->kappa(p[celli], Ts); // kvm
        }
    }

    tkappa.ref().correctBoundaryConditions();

    return tkappa;
}


tmp<volScalarField> liquidFilmThermo::hs
(
    const volScalarField& T,
    const bool withRespectToRef = true
)
const
{
    tmp<volScalarField> ths
    (
        new volScalarField
        (
            IOobject
            (
                type() + ":hs",
                film().time().timeName(),
                film().regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            film().regionMesh(),
            dimensionedScalar(dimEnergy/dimMass, Zero),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    scalarField& hs = ths.ref().primitiveFieldRef();

    if (useReferenceValues_)
    {
        forAll(hs, celli)
        {
            hs[celli] = this->hs(pRef_, Tref_,false);
        }
    }
    else
    {
        const thermoSingleLayer& film = thermoFilm();
        
        const volScalarField& p = film.pPrimary();
        
        forAll(hs, celli)
        {
            const scalar Ts = max(min(T[celli],Tmax_),Tmin_); 
            hs[celli] = this->hs(p[celli], Ts, withRespectToRef); 
        }
    }

    ths.ref().correctBoundaryConditions();

    return ths;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
