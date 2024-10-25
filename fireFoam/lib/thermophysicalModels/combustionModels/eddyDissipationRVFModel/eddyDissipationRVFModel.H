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

Class
    Foam::combustionModels::eddyDissipationRVFModel

Description
    Implementation of: "Modeling of flame extinction/re-ignition in oxygen-reduced 
    environments" by Ren, Zeng, Meredith, Wang, and Dorofeev
    in Proceedings of the Combustion Institute 37(3) 2019.

    Based on "Thermal quenching of mixed eddies in non-premixed flames" 
    by Dorofeev
    in Proceedings of the Combustion Institute 36(2) 2017.

    This is a combustion model based on the eddy dissipaction concept (EDC)
    modified to include local extinction.

SourceFiles
    eddyDissipationRVFModel.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDissipationRVFModel_H
#define eddyDissipationRVFModel_H

#include "singleStepCombustion.H"

#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class eddyDissipationRVFModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class eddyDissipationRVFModel
:
    public singleStepCombustion<ReactionThermo, ThermoType>
{
    // Private data

    //- Model constant
    scalar C_;

    scalar Cd_;

    scalar Cstiff_;

    //- Time to Turnon the Extinction Model
	scalar tExt_;

    scalar Cevap_;

	//- Zeldovich Number at normal-air condition
	scalar ZN_;

	//- Model calibration coefficients
	scalar cKa_;
	scalar cKapa_;

    //- Extra heat loss term for adiabatic flame temperature (default 0)
    scalar XrExt_;

    //- Initial fuel and air stream temperature
	scalar TFuel_;
	scalar TAir_;

    //- Adabatic Flame temperature at normal air condition
	scalar TadAir_;

	//- Flame velocity coefficient
	scalar SLC1_;
	scalar SLC2_;
    //- Effective soot volume fraction at normal air condition
    scalar fvSootAir_;

    //- Critical oxygen volume fraction below which soot vanishes
    scalar O2Soot_;
	//- For adiabatic flame temperature model
	//scalar XrExt_;
	//scalar XO2Air_;
	//scalar rN2_;
	//scalar rCO2_;
	//scalar rH2O_;
    //- Using fixed radiant fraction or variable radiant fraction (default: variable Xr)
    bool fixedXr_;
  
    bool RVFModelActivated_;
    bool variableXr_;


    //- Sergey's Model

    //- Karlovitz Number
    volScalarField Ka_;
    volScalarField KaExt_;
    volScalarField KaMixed_;
    volScalarField KaExtMixed_;

    //- Flame Reactive Volume
    volScalarField RVF_;
    volScalarField PV_;

    //- Flame Extinction Factor (1-RVF_)
    volScalarField FExt_;

    //- Ignition Factor
    volScalarField Fig_;

	//- Energy dissipation rate
	volScalarField epsSgs_;
	volScalarField epsG_;

    //- Laminar Premixed Flame Speed (m/s)
    volScalarField SL_;
    volScalarField SLMixed_;

	//- Reaction rates (W/m3)
	volScalarField dQFuel_;
	volScalarField dQFstar_;

	//- Flame thickness
	volScalarField deltaFlame_;
	volScalarField deltaFlameMixed_;
	//- Heat Loss Term alpha
	volScalarField alphaLoss_;

	//- Expension Ratio
	volScalarField ExpR_;
	volScalarField ExpRMixed_;

	//- Entrainment Ratio
	volScalarField ER_;

	//- Vairable Beta
	volScalarField Beta_;

    //- Gama = 6*G(n)*Kappa0*(Tb/T0)^1.75
	volScalarField Gama_;
	volScalarField GamaMixed_;

    //- Total Extinction Criteria
	volScalarField ExtNumber_;
	volScalarField ExtNumberMixed_;

    //- Adabatic Flame temperature
	volScalarField Tad_;
	volScalarField TadMixed_;

    //- Local oxygen volume fraction
    volScalarField XO2Local_;

    //- Local flame radiant fraction
    volScalarField XrFlame_;

    //- Effective soot volume fraction
    volScalarField fvSoot_;

    //- Radiant power (W/m3): dQ*XrFlame
    volScalarField QdotRad_;
    
    //- Reaction rates (kg/m3/s)
	volScalarField WFstar_;
  
    // Private Member Functions

        //- Return the reciprocal of the turbulent mixing time scale
        tmp<volScalarField> rtTurb() const;

        //- Return the reciprocal of the diffusion time scale
        tmp<volScalarField> rtDiff() const;

        //- No copy construct
        eddyDissipationRVFModel(const eddyDissipationRVFModel&) = delete;

        //- No copy assignment
        void operator=(const eddyDissipationRVFModel&) = delete;


public:

    //- Runtime type information
    TypeName("eddyDissipationRVFModel");


    // Constructors

        //- Construct from components
        eddyDissipationRVFModel
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );


    //- Destructor
    virtual ~eddyDissipationRVFModel();


    // Member Functions

        // Evolution

            //- Correct combustion rate
            virtual void correct();

	    //- Fuel consumption rate matrix                                     
            virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;                       

            //- Sensible enthalpy source term                                        
            //virtual tmp<volScalarField> Sh() const;
            //virtual tmp<volScalarField> ShFuel() const;
            //virtual tmp<volScalarField> ShFstar() const;
            virtual tmp<volScalarField> Qdot() const;
            virtual tmp<volScalarField> QdotFuel() const;
            virtual tmp<volScalarField> QdotFstar() const;

	    void calculateFlameTemperature();
	    void calculateReactiveVolume();

        tmp<volScalarField> rtSpread() const;

        //bool fixedXr();
        //bool modelActivated();
        bool variableXr() const;
        scalar RVFactivationTime() const;


        // I-O

            //- Update properties
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "eddyDissipationRVFModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
