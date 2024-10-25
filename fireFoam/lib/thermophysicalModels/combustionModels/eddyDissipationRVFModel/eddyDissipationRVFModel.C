/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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

#include "eddyDissipationRVFModel.H"
#include "turbulenceModel.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"
#include "fvCFD.H"

#include "thermo.H"
#include "janafThermo.H"
#include "absoluteEnthalpy.H"
#include "perfectGas.H"
#include "etcFiles.H"
#include "IFstream.H"

namespace Foam
{
    typedef species::thermo<janafThermo<perfectGas<specie>>, absoluteEnthalpy>
    thermo;
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
eddyDissipationRVFModel<ReactionThermo, ThermoType>::eddyDissipationRVFModel
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    singleStepCombustion<ReactionThermo, ThermoType>
    (
	    modelType,
	    thermo,
	    turb,
	    combustionProperties
    ),
    C_(readScalar(this->coeffs().lookup("C_EDC"))),
    Cd_(readScalar(this->coeffs().lookup("C_Diff"))),
    Cstiff_(readScalar(this->coeffs().lookup("C_Stiff"))),
    tExt_(this->coeffs().lookupOrDefault("ExtinctionStart", 5.0)),
    Cevap_(this->coeffs().lookupOrDefault("Cevap", 0.5)),

    ZN_(this->coeffs().lookupOrDefault("ZN", 10.0)),
    cKa_(this->coeffs().lookupOrDefault("cKa", 1.0)),
    cKapa_(this->coeffs().lookupOrDefault("cKapa", 1.0)),
    XrExt_(this->coeffs().lookupOrDefault("XrExt", 0.0)),

    TFuel_(this->coeffs().lookupOrDefault("TFuel", 293.15)),
    TAir_(this->coeffs().lookupOrDefault("TAir", 293.15)),
    TadAir_(this->coeffs().lookupOrDefault("TadAir", 2400)),
    SLC1_(this->coeffs().lookupOrDefault("SLC1", 56.0)),
    SLC2_(this->coeffs().lookupOrDefault("SLC2", 11.4)),
    fvSootAir_(this->coeffs().lookupOrDefault("fvSootAir", 2.3)),   //- in PPM
    O2Soot_(this->coeffs().lookupOrDefault("O2Soot", 0.137)),
    fixedXr_(this->coeffs().lookupOrDefault("fixedXr", false)),
    RVFModelActivated_(false),
    variableXr_(false),

    Ka_
    (
        IOobject
        (
            "Ka",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            //IOobject::AUTO_WRITE
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    KaExt_
    (
        IOobject
        (
            "KaExt",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            //IOobject::AUTO_WRITE
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    KaMixed_
    (
        IOobject
        (
            "KaMixed",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    KaExtMixed_
    (
        IOobject
        (
            "KaExtMixed",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    RVF_
    (
        IOobject
        (
            "RVF",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            //IOobject::AUTO_WRITE
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 1.0)
    ),
    PV_
    (
        IOobject
        (
            "PV",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
            //IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    FExt_
    (
        IOobject
        (
            "FExt",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 1.0)
    ),
    Fig_
    (
        IOobject
        (
            "Fig",
            this->mesh().time().timeName(),
            this->mesh(),
            //IOobject::NO_READ,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0),
	    zeroGradientFvPatchScalarField::typeName
    ),
    epsSgs_
    (
        IOobject
        (
            "epsSgs",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimLength*dimLength/dimTime/dimTime/dimTime, 0.0)
    ),
    epsG_
    (
        IOobject
        (
            "epsilon",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimLength*dimLength/dimTime/dimTime/dimTime, 0.0)
    ),
    SL_
    (
        IOobject
        (
            "SL",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            //IOobject::AUTO_WRITE
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimLength/dimTime, 0.0)
    ),
    SLMixed_
    (
        IOobject
        (
            "SLMixed",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimLength/dimTime, 0.0)
    ),
    dQFuel_
    (
        IOobject
        (
            "dQFuel",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            //IOobject::AUTO_WRITE
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    dQFstar_
    (
        IOobject
        (
            "dQFstar",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            //IOobject::AUTO_WRITE
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    deltaFlame_
    (
        IOobject
        (
            "deltaFlame",
            this->mesh().time().timeName(),
            this->mesh(),
            //IOobject::READ_IF_PRESENT,
            IOobject::NO_READ,
            //IOobject::AUTO_WRITE
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    deltaFlameMixed_
    (
        IOobject
        (
            "deltaFlameMixed",
            this->mesh().time().timeName(),
            this->mesh(),
            //IOobject::READ_IF_PRESENT,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    alphaLoss_
    (
        IOobject
        (
            "alphaLoss",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    ExpR_
    (
        IOobject
        (
            "ExpR",
            this->mesh().time().timeName(),
            this->mesh(),
            //IOobject::NO_READ,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 1.0)
    ),
    ExpRMixed_
    (
        IOobject
        (
            "ExpRMixed",
            this->mesh().time().timeName(),
            this->mesh(),
            //IOobject::NO_READ,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 1.0)
    ),
    ER_
    (
        IOobject
        (
            "ER",
            this->mesh().time().timeName(),
            this->mesh(),
            //IOobject::NO_READ,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 1.0)
    ),
    Beta_
    (
        IOobject
        (
            "Beta",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            //IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 10.0)
    ),
    Gama_
    (
        IOobject
        (
            "Gama",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 10.0)
    ),
    GamaMixed_
    (
        IOobject
        (
            "GamaMixed",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 10.0)
    ),
    ExtNumber_
    (
        IOobject
        (
            "ExtNumber",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            //IOobject::AUTO_WRITE
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 10.0)
    ),
    ExtNumberMixed_
    (
        IOobject
        (
            "ExtNumberMixed",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 10.0)
    ),
    Tad_
    (
        IOobject
        (
            "Tad",
            this->mesh().time().timeName(),
            this->mesh(),
            //IOobject::NO_READ,
            IOobject::READ_IF_PRESENT,
            //IOobject::AUTO_WRITE
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimTemperature, 293.0)
    ),
    TadMixed_
    (
        IOobject
        (
            "TadMixed",
            this->mesh().time().timeName(),
            this->mesh(),
            //IOobject::NO_READ,
            IOobject::READ_IF_PRESENT,
            //IOobject::AUTO_WRITE
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimTemperature, 293.0)
    ),
    XO2Local_
    (
        IOobject
        (
            "XO2Local",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            //IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.21)
    ),
    XrFlame_
    (
        IOobject
        (
            "XrFlame",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            //IOobject::READ_IF_PRESENT,
            //IOobject::AUTO_WRITE
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 0.2)
    ),
    fvSoot_
    (
        IOobject
        (
            "fvSoot",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            //IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimless, 2.0)
    ),
    QdotRad_
    (
        IOobject
        (
            "QdotRad",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimEnergy/pow3(dimLength)/dimTime, 0.0)
    ),
    WFstar_
    (
        IOobject
        (
            "WFstar",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    )
{
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
eddyDissipationRVFModel<ReactionThermo, ThermoType>::~eddyDissipationRVFModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationRVFModel<ReactionThermo, ThermoType>::rtTurb() const
{
    return C_*this->turbulence().epsilon()/
              max(this->turbulence().k(),
              dimensionedScalar("SMALL",dimVelocity*dimVelocity,SMALL));
}

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationRVFModel<ReactionThermo, ThermoType>::rtDiff() const
{
    const volScalarField& YO2 = this->thermo().composition().Y("O2");
    const compressible::LESModel& lesModel =
        YO2.db().lookupObject<compressible::LESModel>
	(
	turbulenceModel::propertiesName
	);

    //return Cd_*this->thermo().alpha()/this->rho()/sqr(lesModel.delta());

    dimensionedScalar Df("Df", dimLength, 0.010);
    return Cd_*this->thermo().alpha()/this->rho()/sqr(min(Df,lesModel.delta()));
}

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationRVFModel<ReactionThermo, ThermoType>::rtSpread() const
{
    const volScalarField& YO2 = this->thermo().composition().Y("O2");
    const compressible::LESModel& lesModel =
        YO2.db().lookupObject<compressible::LESModel>
        (
        turbulenceModel::propertiesName
        );

    return SLMixed_/lesModel.delta();
}

template<class ReactionThermo, class ThermoType>
void eddyDissipationRVFModel<ReactionThermo, ThermoType>::correct()
{
    this->wFuel_ ==
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0);

    if (this->active())
    {
        this->singleMixturePtr_->fresCorrect();

        const label fuelI = this->singleMixturePtr_->fuelIndex();

        const volScalarField& YFuel =
            this->thermo().composition().Y()[fuelI];

        const volScalarField& YFstar =
	    this->thermo().composition().Y("Fstar");

        const dimensionedScalar s = this->singleMixturePtr_->s();

        if (this->thermo().composition().contains("O2"))
        {
            const volScalarField& YO2 = this->thermo().composition().Y("O2");

            volScalarField rt(max(rtTurb(),rtDiff()));

	    calculateFlameTemperature();
	    calculateReactiveVolume();

            this->wFuel_ ==
                this->rho()
                * min(YFuel, YO2/s.value())
                / this->mesh().time().deltaT() / Cstiff_
                * (1 - exp(- Cstiff_*this->mesh().time().deltaT() * rt));

            this->WFstar_ ==
                this->rho()
                * min(YFstar, YO2/s.value())
                / this->mesh().time().deltaT() / Cstiff_
                * (1 - exp(- Cstiff_*this->mesh().time().deltaT() * rt));
        }
	    dQFuel_ = this->QdotFuel();
	    dQFstar_ = this->QdotFstar();
        QdotRad_ = XrFlame_*(dQFuel_ + dQFstar_);
    }
}

template<class ReactionThermo, class ThermoType>
bool eddyDissipationRVFModel<ReactionThermo, ThermoType>::read()
{
    if (singleStepCombustion<ReactionThermo, ThermoType>::read())
    {
        this->coeffs().lookup("C") >> C_ ;
        return true;
    }
    else
    {
        return false;
    }
}

template<class ReactionThermo, class ThermoType>
void eddyDissipationRVFModel<ReactionThermo, ThermoType>::calculateReactiveVolume()
{
    const label fuelI = this->singleMixturePtr_->fuelIndex();
    const volScalarField& YFuel = this->thermo().composition().Y()[fuelI];
    const volScalarField& TCellRef = this->thermo().T();
    const volScalarField& YFstar = this->thermo().composition().Y("Fstar");
    const volScalarField& YO2 = this->thermo().composition().Y("O2");

    const volScalarField& ft = this->mesh().template lookupObject<volScalarField>("ft");

    tmp<volScalarField> tnut(this->turbulence().nut());
    const volScalarField& nutRef = tnut();
    epsSgs_ = this->turbulence().epsilon();

    forAll(YO2, cellI)
    {
        //- Initialize parameters for normal fuel
        SL_[cellI] = 0;
        deltaFlame_[cellI] = 0;
        Ka_[cellI] = 0;
        KaExt_[cellI] = 1;
        Gama_[cellI] = 6.0;
        RVF_[cellI] = 0;
        ExtNumber_[cellI] = 0;
        alphaLoss_[cellI] = 0;
        Beta_[cellI] = 20.0;
        FExt_[cellI] = 0;

        //- Initialize parameters for quenched fuel
        SLMixed_[cellI] = 0;
        deltaFlameMixed_[cellI] = 0;
        KaMixed_[cellI] = 0;
        KaExtMixed_[cellI] = 1;
        GamaMixed_[cellI] = 6.0;
        ExtNumberMixed_[cellI] = 0;
        Fig_[cellI] = 0;

        //- Calculate radiative heat loss term
        scalar krad(0.7);
        scalar kai(pow(Tad_[cellI]/TAir_, 1.75)*1.4e-5/0.75);
        epsG_[cellI] = epsSgs_[cellI]*(kai + nutRef[cellI])/(1.0e-8 + nutRef[cellI]);
        scalar BetaFuel(10);
        scalar BetaFstar(10);

        //- Local flame radiant fraction
        XrFlame_[cellI] = 0;
        scalar XrFuel(0);
        scalar XrFstar(0);

        if((YFuel[cellI] > 1.0e-4) && (YO2[cellI] > 1.0e-4))
        {
            SL_[cellI] = max(0.0001,SLC1_*exp(-1000.0*SLC2_/Tad_[cellI]));
            deltaFlame_[cellI] = min(0.1, kai/SL_[cellI]);
            Gama_[cellI] = 6.0*cKapa_*pow((Tad_[cellI]/TAir_), 1.75);

            BetaFuel = max(6.0, (min(20.0,ZN_*sqr(TadAir_/Tad_[cellI])*(Tad_[cellI]-TAir_)/(TadAir_-TAir_))));
            Beta_[cellI] = BetaFuel;
            scalar Xext(((BetaFuel+0.667)+sqrt(sqr(BetaFuel+0.667)-6.667*BetaFuel))/(2*BetaFuel));
            scalar Xaut(((BetaFuel+0.667)-sqrt(sqr(BetaFuel+0.667)-6.667*BetaFuel))/(2*BetaFuel));
            scalar Cka(sqr(ExpR_[cellI]*ZN_)/Gama_[cellI]/exp(BetaFuel));
            scalar Te0(Xext*Tad_[cellI]);

            scalar xo(XO2Local_[cellI]);
            scalar kradgas(233.67*xo*xo*xo-110.62*xo*xo+12.49*xo+0.3849);
            scalar alphaGas(5.33*5.67e-8*kradgas*sqr(deltaFlame_[cellI])*pow(Te0,3)/(cKapa_*kai*(1.2*TAir_/Te0)*530*pow(Te0,0.1131)));

            //-Calculate local radiant fraction
            scalar kradTotal(kradgas + 1.803e-3*fvSoot_[cellI]*Te0);
            alphaLoss_[cellI] = min(1.0,5.33*5.67e-8*kradTotal*sqr(deltaFlame_[cellI])*pow(Te0,3)/(cKapa_*kai*(1.2*TAir_/Te0)*530*pow(Te0,0.1131)));
            XrFuel = min(0.5, alphaLoss_[cellI]/(1.0 + alphaLoss_[cellI]));

            Ka_[cellI] = min(10., cKa_*sqr(deltaFlame_[cellI])*sqrt(epsG_[cellI])/pow(kai,1.5));
            KaExt_[cellI] = pow(max(0.0,Cka*(1.0-Xext)*exp(BetaFuel*Xext)/pow(Xext,5.0/3.0)-alphaLoss_[cellI]*pow(Xext,4.0/3.0)), 1.5);
            ExtNumber_[cellI] = KaExt_[cellI] - 1.0/Xext;

            if(ExtNumber_[cellI] < 0)	//-Purely thermal quenching
            {
                RVF_[cellI] = 0;
                FExt_[cellI] = 1;
            }
            else
            {
                if(Ka_[cellI] > KaExt_[cellI])	//-Extinction due to high strain and thermal
                {
                    RVF_[cellI] = 0;
                    FExt_[cellI] = 1;
                }
                else if(Ka_[cellI] < 1.1)	//-Healthy flame
                {
                    RVF_[cellI] = 1.0;
                    FExt_[cellI] = 0.0;
                }
                else
                {
                    if(KaExt_[cellI] - Ka_[cellI] < 0.05)   //-Close to quenching limit
                    {
                   	    RVF_[cellI] = 0.3;
                        FExt_[cellI] = 0.7;
                    }
                    else			//-Partial extinction
                    {
                        scalar X1(Xext-0.05);
                        scalar X2(Xext+0.05);
                        scalar X3(Xext-0.05);
                        scalar iterKa(0);
                        while(iterKa < 10)
                        {
                            iterKa = iterKa + 1;
                            X1 = log(X1*(pow(Ka_[cellI]*X1,0.667)+alphaLoss_[cellI]*sqr(X1))/(Cka*(1-X1)))/BetaFuel;
                            scalar YX2(exp(BetaFuel*X2)*Cka*pow(X2,0.333)/(pow(Ka_[cellI],0.667)+alphaLoss_[cellI]*pow(X2,1.333)));
                            X2 = (sqrt(sqr(YX2)+4.0*YX2)-YX2)/2.0;
                            X3 = log(X3*(1+alphaLoss_[cellI]*sqr(X3))/(Cka*(1-X3)))/BetaFuel;
                        }
                        //RVF_[cellI] = (pow(X3/X1, 5.0) - pow(X3/X2, 5.0))/(1.0 - pow(X3, 5.0));
                        RVF_[cellI] = max(0.0, min(1.0, (pow(X3/X1, 5.0) - pow(X3, 5.0))/(1.0 - pow(X3, 5.0))));
                        FExt_[cellI] = 1.0 - RVF_[cellI];
                    }
                }
            }

            if(YO2[cellI] < 1.0e-4)	//-Let pure fuel be flammable
            {
                RVF_[cellI] = 1.0;
                FExt_[cellI] = 0.0;
            }
        }

        if(YFstar[cellI] > 1.0e-4)
        {
            scalar kaiMixed(pow(TadMixed_[cellI]/TAir_, 1.75)*1.4e-5/0.75);
            BetaFstar = max(6.0, (min(20.0,ZN_*sqr(TadAir_/TadMixed_[cellI])*(TadMixed_[cellI]-TAir_)/(TadAir_-TAir_))));
            SLMixed_[cellI] = max(0.0001,SLC1_*exp(-1000.0*SLC2_/TadMixed_[cellI]));
            deltaFlameMixed_[cellI] = min(0.1, kaiMixed/SLMixed_[cellI]);
            GamaMixed_[cellI] = 6.0*cKapa_*pow((TadMixed_[cellI]/TAir_), 1.75);

            scalar XextMixed(((BetaFstar+0.667)+sqrt(sqr(BetaFstar+0.667)-6.667*BetaFstar))/(2*BetaFstar));
            scalar XautMixed(((BetaFstar+0.667)-sqrt(sqr(BetaFstar+0.667)-6.667*BetaFstar))/(2*BetaFstar));
            scalar CkaMixed(sqr(ExpRMixed_[cellI]*ZN_)/GamaMixed_[cellI]/exp(BetaFstar));
            scalar Te0Mixed(XextMixed*TadMixed_[cellI]);
            scalar CalphaMixed(5.33*5.67e-8*krad*sqr(deltaFlameMixed_[cellI])*pow(Te0Mixed,3)/(cKapa_*kaiMixed*(1.2*TAir_/Te0Mixed)*530*pow(Te0Mixed,0.1131)));

            KaMixed_[cellI] = min(10., cKa_*sqr(deltaFlameMixed_[cellI])*sqrt(epsG_[cellI])/pow(kaiMixed,1.5));
            KaExtMixed_[cellI] = pow(max(0.0,CkaMixed*(1.0-XextMixed)*exp(BetaFstar*XextMixed)/pow(XextMixed,5.0/3.0)-CalphaMixed*pow(XextMixed,4.0/3.0)) , 1.5);
            ExtNumberMixed_[cellI] = KaExtMixed_[cellI] - 1.0/XextMixed;

            //- Simplified treatment for re-ignition
            if(ExtNumberMixed_[cellI] < 0)
            {
        	    Fig_[cellI] = 0;
            }
            else
            {
                if(KaMixed_[cellI] > KaExtMixed_[cellI])
                {
                    Fig_[cellI] = 0;
                }
                else
                {
                    Fig_[cellI] = 1.0;
                }
            }
        }

        //- Local flame radiant fraction contributed from normal and quenched fuel
        XrFlame_[cellI] = (XrFuel*YFuel[cellI]*RVF_[cellI] + XrFstar*YFstar[cellI]*Fig_[cellI])/(YFuel[cellI]*RVF_[cellI] + YFstar[cellI]*Fig_[cellI] + 1.0e-8);

        //if(this->mesh().time().value() < tExt_)	//-Activation flame extinction model after tExt
        if(!RVFModelActivated_)
        {
            RVF_[cellI] = 1.0;
            Fig_[cellI] = 1.0;
            //XrFlame_[cellI] = 0.22;
        }
    }

    if(this->mesh().time().value() < tExt_)	//-Activation flame extinction model after tExt
    {
        RVFModelActivated_=false;
    }
    else
    {
        RVFModelActivated_=true;
    }

    //if(fixedXr_)
    //{
    //    Info<<"Using Fixed Radiant Fraction!"<<endl;
    //    forAll(XrFlame_,cellI)
    //    {
    //        //- Hard coded for bio-fuel
    //        XrFlame_[cellI] = 0.225;
    //    }
    //}
    if(RVFModelActivated_ && !fixedXr_)
    {
        Info<<"Using RVF soot emission model!"<<endl;
        variableXr_=true;
    }
    else
    {
        Info<<"Using constant emission model!"<<endl;
        variableXr_=false;
    }
}

template<class ReactionThermo, class ThermoType>
void eddyDissipationRVFModel<ReactionThermo, ThermoType>::calculateFlameTemperature()
{
    //- Get species mass fraction
    const label fuelI = this->singleMixturePtr_->fuelIndex();
    const volScalarField& YFuel = this->thermo().composition().Y()[fuelI];
    const volScalarField& YO2 = this->thermo().composition().Y("O2");
    const volScalarField& YN2 = this->thermo().composition().Y("N2");
    const volScalarField& YCO2 = this->thermo().composition().Y("CO2");
    const volScalarField& YH2O = this->thermo().composition().Y("H2O");
    const volScalarField& YFstar = this->thermo().composition().Y("Fstar");

    const volScalarField& TCellRef = this->thermo().T();
    const volScalarField& pCellRef = this->thermo().p();
    const dimensionedScalar qF(this->singleMixturePtr_->qFuel());
    const dimensionedScalar s = this->singleMixturePtr_->s();
    const volScalarField rhoCell(this->thermo().rho());


    fileName thermoDataFileName("./constant/thermo.compressibleGas");
    IFstream thermoDataFile(thermoDataFileName);
    dictionary thermoData(thermoDataFile);
    thermo TFuel
    (
        "TFuel",
	    thermo(thermoData.subDict(this->thermo().composition().species()[fuelI]))
    );
    thermo TO2
    (
        "TO2",
	    thermo(thermoData.subDict("O2"))
    );
    thermo TN2
    (
        "TN2",
	    thermo(thermoData.subDict("N2"))
    );
    thermo TCO2
    (
        "TCO2",
	    thermo(thermoData.subDict("CO2"))
    );
    thermo TH2O
    (
        "TH2O",
	    thermo(thermoData.subDict("H2O"))
    );

    //- Get Mspecies/Mfuel from reaction equation
    scalar rCO2(this->singleMixturePtr_->specieStoichCoeffs()
            [this->thermo().composition().species()["CO2"]]);
    scalar rH2O(this->singleMixturePtr_->specieStoichCoeffs()
            [this->thermo().composition().species()["H2O"]]);
    scalar rN2(this->singleMixturePtr_->specieStoichCoeffs()
            [this->thermo().composition().species()["N2"]]);

    //- Get Spray Info
    //const volScalarField sprayDensity =
    //   	this->mesh().template lookupObject<volScalarField>("rhoSpray");
    volScalarField sprayDensity(this->mesh().template lookupObject<volScalarField>("rhoSpray"));
    if(this->mesh().template foundObject<volScalarField>("rhoSprayMean_MA"))
    {
        Info<<"Using moving averaged spray density."<<endl;
        sprayDensity = this->mesh().template lookupObject<volScalarField>("rhoSprayMean_MA");
    }

    PV_ = (YCO2*(1.0+rH2O/rCO2)+SMALL)/(YCO2*(1.0+rH2O/rCO2)+SMALL + min(YFuel,YO2/s.value())*(1.0+s.value()));

    forAll(YFuel, cellI)
    {
	    ER_[cellI] = 0.0;
	    ExpR_[cellI] = 1.0;
	    ExpRMixed_[cellI] = 1.0;
	    scalar pValue(pCellRef[cellI]);
	    scalar TValue(TCellRef[cellI]);

	    //- Calculate local original O2 mole fraction (convert combustion product to original species)
        scalar O2Total(YCO2[cellI]/TCO2.W()+0.5*YCO2[cellI]*(rH2O/rCO2)/TH2O.W()+YO2[cellI]/TO2.W());
        XO2Local_[cellI]=min(0.25, max(O2Total/(O2Total + YN2[cellI]/TN2.W() + 1.0e-6),1.0e-6));
        scalar XO2Air(XO2Local_[cellI]);
        rN2  = s.value()*(1-XO2Air)*TN2.W()/XO2Air/TO2.W();
	    scalar Yspray(max(min(Cevap_*sprayDensity[cellI]/rhoCell[cellI],1.0),0.0));

	    //- Calculate adiabatic flame temperature for normal fuel
	    if((YFuel[cellI] > 1.0e-4) && (YO2[cellI] > 1.0e-4))
	    {
	        scalar YN2R(min(YN2[cellI], YO2[cellI]*TN2.W()*(1.0 - XO2Air)/(TO2.W()*XO2Air)));
	        scalar YN2P(max(0.0,YN2[cellI]-YN2R));
	        scalar YFuelR(min(YFuel[cellI],YO2[cellI]/s.value()));
	        scalar MFO(TFuel.W()/TO2.W());
	        scalar MFN(TFuel.W()/TN2.W());
	        ER_[cellI] = YFuelR*(1.0+s.value()*MFO+s.value()*MFN*YN2R/YO2[cellI])
	                     /
	    		 (YFuel[cellI]+MFO*YO2[cellI]+MFN*YN2[cellI]);
	        scalar ERValue(ER_[cellI]);
	        scalar CoffN2(ERValue*YN2P+rN2*YFuelR);
	        scalar CoffCO2(ERValue*YCO2[cellI]+rCO2*YFuelR);
	        scalar CoffH2O(ERValue*(YH2O[cellI]+Yspray)+rH2O*YFuelR);
	        scalar CoffFstar(ERValue*YFstar[cellI]);
	        scalar CoffO2(0);
	        scalar RHS1(YFuelR*qF.value()*(1.0-XrExt_));
	        scalar RHS2(min(RHS1,ERValue*Yspray*3.0e6));
	        scalar RHS3(ERValue*
	    		          (YO2[cellI]*TO2.Hs(pValue,TValue)
	    	              +YN2[cellI]*TN2.Hs(pValue,TValue)
	    	              +YCO2[cellI]*TCO2.Hs(pValue,TValue)
	    	              +YH2O[cellI]*TH2O.Hs(pValue,TValue)
	    	              +(YFuel[cellI]+YFstar[cellI])*TFuel.Hs(pValue,TValue)
	    	              )
	                   );
	        scalar RHS4(
	    	              (YFuelR - ERValue*YFuel[cellI])*TFuel.Hs(pValue,TFuel_)
	    	              + (s.value()*YFuelR - ERValue*YO2[cellI])*TO2.Hs(pValue,TAir_)
	    	              + (rN2*YFuelR - ERValue*YN2R)*TN2.Hs(pValue,TAir_)
	                   );

	        scalar RHSall(RHS1 - RHS2 + RHS3 + RHS4);

	        scalar Test(TValue);
	        scalar Tnew(TValue);
	        int iter(0);
	        do
	        {
	            Test = Tnew;
	            scalar CpEff(
	                         CoffN2*TN2.Cp(pValue,Test)
	            	        +CoffCO2*TCO2.Cp(pValue,Test)
	            	        +CoffH2O*TH2O.Cp(pValue,Test)
	            	        +CoffFstar*TFuel.Cp(pValue,Test)
	            	        +CoffO2*TO2.Cp(pValue,Test)
	                        );
	            scalar LHSest(
	                         CoffN2*TN2.Hs(pValue,Test)
	            	        +CoffCO2*TCO2.Hs(pValue,Test)
	            	        +CoffH2O*TH2O.Hs(pValue,Test)
	            	        +CoffFstar*TFuel.Hs(pValue,Test)
	            	        +CoffO2*TO2.Hs(pValue,Test)
	                        );
	            Tnew = Test + (RHSall - LHSest)/CpEff;
	            if (iter++ > 10)
	            {
	                FatalErrorInFunction
	                    << "Maximum number of iteration exceeded in Tad model" << Tnew <<endl
	                    <<"Local XO2: "<<XO2Air<<tab<<cellI<<tab<<iter<<endl
	                    <<"Coeffs: "<<CoffN2<<tab<<CoffH2O<<tab<<CoffCO2<<tab<<CoffO2<<endl
	                    <<"RHSall: "<<RHSall<<tab<<LHSest<<tab<<CpEff<<endl
	                    << abort(FatalError);
	            }
	            if ((Tnew > 5000) || (Tnew < 200))
	            {
	                FatalErrorInFunction
	                    <<"Tad exceed range of [200, 5000]: " << Tnew <<endl
	                    <<"Local XO2: "<<XO2Air<<tab<<cellI<<tab<<iter<<endl
	                    <<"O2 (total 1) & N2: "<<O2Total<<tab<<YN2[cellI]<<endl
	                    <<"Coeffs: "<<CoffN2<<tab<<CoffH2O<<tab<<CoffCO2<<tab<<CoffO2<<endl
	           	        <<"RHSall: "<<RHSall<<tab<<LHSest<<tab<<CpEff<<endl
	           	        <<"Yspray: "<<Yspray<<endl
	           	        <<"ERValue: "<<ERValue<<endl
	                    << abort(FatalError);
	            }
	        } while (mag(Tnew - Test) > 5.0);

	        Tad_[cellI] = Tnew;
	        ExpR_[cellI] = (Tnew/TAir_)
	                      *( (rCO2/TCO2.W()+rH2O/TH2O.W()+rN2/TN2.W()) /
	                         (1.0/TFuel.W()+s.value()/TO2.W()+rN2/TN2.W())
	                       );
	    }
	    else
	    {
	        Tad_[cellI] = TValue;
	    }

	    //- Calculate adiabatic flame temperature for mixed fuel
	    if((YFstar[cellI] > 1.0e-4) && (YO2[cellI] > 1.0e-4))
	    {
	        scalar YFuelR(min(YFstar[cellI],YO2[cellI]/s.value()));
	        scalar CoffN2(YN2[cellI]);
	        scalar CoffH2O(YH2O[cellI]+rH2O*YFuelR+Yspray);
	        scalar CoffCO2(YCO2[cellI]+rCO2*YFuelR);
	        scalar CoffFstar(YFstar[cellI] - YFuelR);
	        scalar CoffO2(YO2[cellI] - s.value()*YFuelR);

	        scalar RHS1(YFuelR*qF.value()*(1.0-XrExt_));
	        scalar RHS2(min(RHS1,Yspray*3.0e6));
	        scalar RHS3( YO2[cellI]*TO2.Hs(pValue,TValue)
	                   + YN2[cellI]*TN2.Hs(pValue,TValue)
	                   + YCO2[cellI]*TCO2.Hs(pValue,TValue)
	                   + YH2O[cellI]*TH2O.Hs(pValue,TValue)
	                   + YFstar[cellI]*TFuel.Hs(pValue,TValue));
	        scalar RHS4(0);

	        scalar RHSall(RHS1 - RHS2 + RHS3 + RHS4);

	        scalar Test(TValue);
	        scalar Tnew(TValue);
	        int iter(0);
	        do
	        {
	    	    Test = Tnew;
	    	    scalar CpEff(
	    	                 CoffN2*TN2.Cp(pValue,Test)
	    	                +CoffCO2*TCO2.Cp(pValue,Test)
	    	                +CoffH2O*TH2O.Cp(pValue,Test)
	    	                +CoffFstar*TFuel.Cp(pValue,Test)
	    	                +CoffO2*TO2.Cp(pValue,Test)
	    	                );
	    	    scalar LHSest(
	    	                 CoffN2*TN2.Hs(pValue,Test)
	    	                +CoffCO2*TCO2.Hs(pValue,Test)
	    	                +CoffH2O*TH2O.Hs(pValue,Test)
	    	                +CoffFstar*TFuel.Hs(pValue,Test)
	    	                +CoffO2*TO2.Hs(pValue,Test)
	    	                );
	    	    Tnew = Test + (RHSall - LHSest)/CpEff;
	    	    if (iter++ > 10)
	    	    {
	    	        FatalErrorInFunction
	    	            << "Maximum number of iteration exceeded in Tad model"
	    	   	        << abort(FatalError);
	    	    }
	    	    if ((Tnew > 5000) || (Tnew < 200))
	    	    {
	    	        FatalErrorInFunction
	    	            <<"Tad exceed range of [200, 5000]: " << Tnew <<endl
	    	            <<"Local XO2: "<<XO2Air<<tab<<cellI<<tab<<iter<<endl
	    	            <<"O2 (total 2) & N2: "<<O2Total<<tab<<YN2[cellI]<<endl
	    	            <<"Coeffs: "<<CoffN2<<tab<<CoffH2O<<tab<<CoffCO2<<tab<<CoffO2<<endl
	    	            <<"RHSall: "<<RHSall<<tab<<LHSest<<tab<<CpEff<<endl
	           	        <<"Yspray: "<<Yspray<<endl
	    	            << abort(FatalError);
	    	    }
	        } while (mag(Tnew - Test) > 5.0);

	        TadMixed_[cellI] = Tnew;
	        ExpRMixed_[cellI] = (Tnew/TAir_)
	    	                   *( (rCO2/TCO2.W()+rH2O/TH2O.W()+rN2/TN2.W()) /
	    	                      (1.0/TFuel.W()+s.value()/TO2.W()+rN2/TN2.W())
	    	                    );
	    }
	    else
	    {
	        TadMixed_[cellI] = TValue;
	    }

        //- Local soot volume fraction near flame sheet
        scalar XO2Eff(0);
        if(YFuel[cellI] > 1.0e-4)
        {
            XO2Eff = max(0, 0.21 - (TadAir_ - Tad_[cellI])/8055);
        }
        else if(YFstar[cellI] > 1.0e4)
        {
            XO2Eff = max(0, 0.21 - (TadAir_ - TadMixed_[cellI])/8055);
        }
        else
        {
            XO2Eff = 0;
        }
        fvSoot_[cellI] = max(0, fvSootAir_*(XO2Eff - O2Soot_)/(0.21-O2Soot_));
    }
}

template<class ReactionThermo, class ThermoType>
tmp<fvScalarMatrix>
eddyDissipationRVFModel<ReactionThermo, ThermoType>::R
(
    volScalarField& Y
) const
{
    const label specieI = this->thermo().composition().species()[Y.name()];
    const label fuelI = this->singleMixturePtr_->fuelIndex();
    if(specieI == fuelI)
    {
        volScalarField wSpecie
        (
            this->wFuel_*this->singleMixturePtr_->specieStoichCoeffs()[specieI]
        );
        return wSpecie + fvm::Sp(0.0*wSpecie, Y);
    }
    else if(Y.name() == "Fstar")
    {
        volScalarField wSpecie
        (
           (1-RVF_)*this->wFuel_ - Fig_*this->WFstar_
        );
        return wSpecie + fvm::Sp(0.0*wSpecie, Y);
    }
    else
    {
        volScalarField wSpecie
        (
           (RVF_*this->wFuel_ + Fig_*this->WFstar_)
	       *this->singleMixturePtr_->specieStoichCoeffs()[specieI]
        );
        return wSpecie + fvm::Sp(0.0*wSpecie, Y);
    }
}


template<class ReactionThermo, class ThermoType>
tmp<volScalarField>
eddyDissipationRVFModel<ReactionThermo, ThermoType>::Qdot() const
{
    const label fuelI = this->singleMixturePtr_->fuelIndex();
    volScalarField& YFuel =
        const_cast<volScalarField&>(this->thermo().composition().Y(fuelI));

    const label indexFstar(this->thermo().composition().species()["Fstar"]);
    volScalarField& YFstar =
        const_cast<volScalarField&>(this->thermo().composition().Y(indexFstar));

    return -this->singleMixturePtr_->qFuel()*((R(YFuel) & YFuel) + (R(YFstar) & YFstar));
}

template<class ReactionThermo, class ThermoType>
tmp<volScalarField>
eddyDissipationRVFModel<ReactionThermo, ThermoType>::QdotFuel() const
{
    const label fuelI = this->singleMixturePtr_->fuelIndex();
    volScalarField& YFuel =
        const_cast<volScalarField&>(this->thermo().composition().Y(fuelI));

    return -this->singleMixturePtr_->qFuel()*(R(YFuel) & YFuel)*RVF_;
}

template<class ReactionThermo, class ThermoType>
tmp<volScalarField>
eddyDissipationRVFModel<ReactionThermo, ThermoType>::QdotFstar() const
{
    const label fuelI = this->singleMixturePtr_->fuelIndex();
    volScalarField& YFuel =
        const_cast<volScalarField&>(this->thermo().composition().Y(fuelI));

    const label indexFstar(this->thermo().composition().species()["Fstar"]);
    volScalarField& YFstar =
        const_cast<volScalarField&>(this->thermo().composition().Y(indexFstar));

    return -this->singleMixturePtr_->qFuel()*
                  ((R(YFuel) & YFuel)*(1-RVF_)+(R(YFstar) & YFstar));
}

//template<class ReactionThermo, class ThermoType>
//bool eddyDissipationRVFModel<ReactionThermo, ThermoType>::fixedXr()
//{
//    return fixedXr_;
//}
//template<class ReactionThermo, class ThermoType>
//bool eddyDissipationRVFModel<ReactionThermo, ThermoType>::modelActivated()
//{
//    return RVFModelActivated_;
//}
template<class ReactionThermo, class ThermoType>
bool eddyDissipationRVFModel<ReactionThermo, ThermoType>::variableXr() const
{
    return variableXr_;
}
template<class ReactionThermo, class ThermoType>
scalar eddyDissipationRVFModel<ReactionThermo, ThermoType>::RVFactivationTime() const
{
    return tExt_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
