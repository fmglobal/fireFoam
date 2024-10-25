Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

#include "readAdditionalThermo.H" // kvm

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField& p = thermo.p();

Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Moved from phrgh.H
volScalarField ph_rgh 
(
    IOobject
    (
        "ph_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

#include "createMRF.H"


Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating combustion model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> combustion
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
#include "readpRef.H"

volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

thermo.correct(); // kvm
rho = thermo.rho(); // kvm

mesh.setFluxRequired(p_rgh.name());

#include "phrghEqn.H"


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
);

#include "createDpdt.H"

#include "createK.H"


volScalarField kappa 
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.Cp()*thermo.alpha()
);

volScalarField cp 
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.Cp()
);

singleStepReactingMixture<gasHThermoPhysics>& singleMixture 
(
    dynamic_cast<singleStepReactingMixture<gasHThermoPhysics>&>
    (thermo)
);

// stoichiometric O2 to fuel ratio
scalar s(singleMixture.s().value()); 

// stoichiometric air to fuel ratio
scalar stoicRatio(singleMixture.stoicRatio().value()); 

// heat of combustion [J/kg]
scalar qFuel(singleMixture.qFuel().value());  

label fuelIndex(singleMixture.fuelIndex()); 
const label inertIndex(composition.species()[inertSpecie]);

const volScalarField* O2Ptr = &Y[inertIndex]; //dummy ptr when O2 is not available 
if (thermo.composition().contains("O2"))  
{
    O2Ptr = &thermo.composition().Y("O2");
}
const volScalarField& O2 = *O2Ptr;

const volScalarField* CO2Ptr = &Y[inertIndex]; //dummy ptr when O2 is not available 
if (thermo.composition().contains("CO2")) 
{
    CO2Ptr = &thermo.composition().Y("CO2");
}
const volScalarField& CO2 = *CO2Ptr;

const volScalarField& fu = Y[fuelIndex]; 

scalar YO2Inf = 0.23301; //hardcode for now

// Fuel molecular weight
scalar Wu_ = singleMixture.speciesData()[fuelIndex].W(); 
// Fuel enthalpy of formation
scalar Hu_ = singleMixture.speciesData()[fuelIndex].hc(); 

// compute stoichiometric mixture fraction
scalar ftSt = 1.0 / ( 1.0 + stoicRatio );
Info << "stoichiometric mixture fraction is = " << ftSt << nl; 

// For vizualizing the smoke, zero val BC needed 
volScalarField smoke
(
    IOobject
    (
        "smoke",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (CO2),
    fixedValueFvPatchScalarField::typeName
);



// create fileds for mixture fraction
volScalarField ft 
(
    IOobject
    (
        "ft",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (fu*s-O2+YO2Inf)/(s+YO2Inf)
);


#include "createClouds.H"
#include "createSurfaceFilmModel.H"
#include "createPyrolysisModel.H"
#include "createRadiationModel.H"
#include "createFvOptions.H"
