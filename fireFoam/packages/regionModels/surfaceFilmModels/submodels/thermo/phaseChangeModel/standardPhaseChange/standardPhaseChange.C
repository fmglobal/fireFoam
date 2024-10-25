/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "standardPhaseChange.H"
#include "addToRunTimeSelectionTable.H"
#include "thermoSingleLayer.H"
#include "zeroField.H"
#include "specie.H" // kvm
#include "heatTransferModel.H" // kvm
#include "filmRadiationModel.H" // kvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardPhaseChange, 0);

addToRunTimeSelectionTable
(
    phaseChangeModel,
    standardPhaseChange,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar standardPhaseChange::Sh
(
    const scalar Re,
    const scalar Sc
) const
{

    if (Re < 5.0E+05)
    {
        return 0.664*sqrt(Re)*cbrt(Sc);
    }
    else
    {
        return 0.037*pow(Re, 0.8)*cbrt(Sc);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardPhaseChange::standardPhaseChange
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    phaseChangeModel(typeName, film, dict),
    deltaMin_(coeffDict_.get<scalar>("deltaMin")),
    L_(coeffDict_.get<scalar>("L")),
    TbFactor_(coeffDict_.getOrDefault<scalar>("TbFactor", 1.1)),
    YInfZero_(coeffDict_.getOrDefault<Switch>("YInfZero", false)),
    scaling_(coeffDict_.getOrDefault<scalar>("scaling",1.0)), // kvm
    dMassFluxMax_(coeffDict_.getOrDefault<scalar>("dMassFluxMax",1)), // ak - limit by mass flux kg/m2/s
    dMassMax_(Zero)                                                      // , to limit mesh-dependency
{ // kvm
    Info << "Mass transfer convective scaling set to " << scaling_ << nl; // kvm


    // AK set dMassMax from dMassFluxMax_ and face areas

    const fvMesh& mesh = film.regionMesh();
    const surfaceScalarField& magSf = mesh.magSf();

    scalar filmArea(Zero);

    for (label i=0; i<film.intCoupledPatchIDs().size(); i++)
    {
        label patchi = film.intCoupledPatchIDs()[i];

        filmArea += gSum(magSf.boundaryField()[patchi]);
    }

    scalar  nCells = returnReduce(mesh.cells().size(), sumOp<label>());

    scalar averageFaceSize = filmArea/nCells;

    dMassMax_ = dMassFluxMax_*averageFaceSize;

//    Info << "Total film area: " << filmArea << endl;
//    Info << "Average face area on film: " << averageFaceSize << endl;
//    Info << "Total number of cells in film: " << nCells << endl;
//    Info << "dMassMax: " << dMassMax_ << endl;

} 



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class YInfType>
void standardPhaseChange::correctModel
(
    const scalar dt,
    scalarField& availableMass,
    scalarField& dMass,
    scalarField& dEnergy,
    YInfType YInf
)
{
    const thermoSingleLayer& film = filmType<thermoSingleLayer>();

    // Set local thermo properties
    const SLGThermo& thermo = film.thermo();
    const filmThermoModel& filmThermo = film.filmThermo();
    const label vapId = thermo.carrierId(filmThermo.name());

    // Retrieve fields from film model
    const scalarField& delta = film.delta();
    const scalarField& pInf = film.pPrimary();
    const scalarField& T = film.T();
    const scalarField& kappa = film.kappa(); // kvm
    const scalarField& rho = film.rho();
    const scalarField& TInf = film.TPrimary(); // kvm
    const scalarField& rhoInf = film.rhoPrimary();
    const scalarField& muInf = film.muPrimary();
    const scalarField& magSf = film.magSf();
    const scalarField hInf(film.htcs().h()); // kvm
    const vectorField dU(film.UPrimary() - film.Us());
    const filmRadiationModel& radiation = film.radiation(); // kvm
    const scalarField  Shs(radiation.Shs());
    const scalarField limMass
    (
        max(scalar(0), availableMass - deltaMin_*rho*magSf)
    );
    
    static scalarField dMassPrev(dMass.size(),0.0); // kvm

    scalar boilingMass = 0.0; // kvm
    scalar evaporationMass = 0.0; // kvm

    forAll(dMass, celli)
    {

        if (delta[celli] > deltaMin_)
        {
            // Cell pressure [Pa]
            const scalar pc = max(min(pInf[celli],101325*1.1),101325*0.9); // kvm

            // Calculate the boiling temperature
            const scalar Tb = filmThermo.Tb(pc);
            // const scalar Tb = 374.8;

            // estimate surface temperature based on energy balance
            const scalar d2k = delta[celli]/2.0/kappa[celli]; // kvm
            const scalar Tf = T[celli]; // kvm
            const scalar Tsurface =  // kvm
                (Tf+d2k*(hInf[celli]*TInf[celli]+Shs[celli])) // kvm
                /(1+d2k*hInf[celli]); // kvm

            // Local temperature - impose lower limit of 200 K for stability
            const scalar Tloc = min(TbFactor_*Tb, max(200.0, Tsurface)); // kvm

            // Saturation pressure [Pa]
            const scalar pSat = filmThermo.pv(pc, Tloc);

            // Latent heat [J/kg]
            const scalar hVap = filmThermo.hl(pc, Tloc);

            // Calculate mass transfer
            if (pSat >= 0.95*pc)
            {
                // Boiling
                const scalar Cp = filmThermo.Cp(pc, Tloc);
                // const scalar Tcorr = 0.25*max(0.0, Tsurface - Tb);
                const scalar Tcorr = 0.25*max(0.0, Tloc - Tb); // kvm
                const scalar qCorr = limMass[celli]*Cp*(Tcorr);
                dMass[celli] = qCorr/hVap;
                boilingMass += dMass[celli]; // kvm
            }
            else
            {
                // Primary region density [kg/m3]
                const scalar rhoInfc = rhoInf[celli];

                // Primary region viscosity [Pa.s]
                const scalar muInfc = muInf[celli];

                // Reynolds number
                scalar Re = rhoInfc*mag(dU[celli])*L_/muInfc; // kvm
                // limit Re to reasonable values, prevent runaway vaporization
                Re = min(1e6,Re); // kvm

                // molecular weight of vapour [kg/kmol]
                const scalar Wvap = thermo.carrier().W(vapId);

                // molecular weight of liquid [kg/kmol]
                const scalar Wliq = filmThermo.W();

                // Vapour mass fraction at interface
                const scalar Ys = Wliq*pSat/(Wliq*pSat + Wvap*(pc - pSat));

                // Vapour diffusivity [m2/s]
                const scalar Dab = filmThermo.D(pc, Tloc);

                // Schmidt number
                scalar Sc = muInfc/(rhoInfc*(Dab + ROOTVSMALL)); // kvm
                if (Sc < 0.01)         // kvm     
                {                      // kvm
                    Sc = max(Sc,0.01); // kvm
                }                      // kvm

                // Sherwood number
                const scalar Sh = this->Sh(Re, Sc);

                // mass transfer coefficient [m/s]
                const scalar hm = scaling_*Sh*Dab/(L_ + ROOTVSMALL); // kvm
                // const scalar hm = scaling_;

                // add mass contribution to source
                dMass[celli] =
                    dt*magSf[celli]*rhoInfc*hm*(Ys - YInf[celli])/(1.0 - Ys);
                evaporationMass += dMass[celli]; // kvm
            }

            dMass[celli] = min(limMass[celli], max(0.0, dMass[celli]));
            // AK: commenting out under-relaxation. It messes things up when timestep is not small
            // try to under-relax the vaporization rate
            //dMass[celli] = min(dMassPrev[celli]*1.05+SMALL, max(0.0, dMass[celli])); // kvm
            // Evaporation limited by max flux value 
            dMass[celli] = min(dMassMax_, dMass[celli]);
    
            dEnergy[celli] = dMass[celli]*hVap;
            dMassPrev[celli] = dMass[celli]; // kvm
        }
    }
    reduce(boilingMass,sumOp<scalar>()); // kvm
    reduce(evaporationMass,sumOp<scalar>()); // kvm
    Info << "boiling fraction " << boilingMass/(boilingMass+evaporationMass+SMALL) << nl; // kvm
}


void standardPhaseChange::correctModel
(
    const scalar dt,
    scalarField& availableMass,
    scalarField& dMass,
    scalarField& dEnergy,
    const label layer  // Alex
)
{
    if (YInfZero_)
    {
        correctModel(dt, availableMass, dMass, dEnergy, zeroField());
    }
    else
    {
        const thermoSingleLayer& film = filmType<thermoSingleLayer>();
        const label vapId = film.thermo().carrierId(film.filmThermo().name());
        const scalarField& YInf = film.YPrimary()[vapId];

        correctModel(dt, availableMass, dMass, dEnergy, YInf);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
