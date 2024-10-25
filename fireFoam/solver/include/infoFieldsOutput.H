
    // create fileds for patch integration of HRR
    surfaceScalarField HRR_fu
    ( 
        IOobject
        (
            "HRR_fu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi*fvc::interpolate(fu)*qFuel
    );

    // for check ft conservation
    surfaceScalarField phiFt
    ( 
        IOobject
        (
            "phiFt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi
    );

    surfaceScalarField phiO2
    ( 
        IOobject
        (
            "phiO2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi
    );

    surfaceScalarField phiCO2
    ( 
        IOobject
        (
            "phiCO2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi
    );

    // create fileds for surface integration of sensible enthalpy flux
    surfaceScalarField phiHs
    (
        IOobject
        (
            "phiHs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi*fvc::interpolate(thermo.he()) 
    );

    // create fileds for surface integration of chemical enthalpy flux
    volScalarField hc(thermo.hc()); 
    surfaceScalarField phiHc
    (
        IOobject
        (
           "phiHc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiHs
    );

    // create fileds for surface integration of total enthalpy flux
    volScalarField h(hc + thermo.he());
    surfaceScalarField phiH
    (
        IOobject
        (
            "phiH",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiHs
    );
    surfaceScalarField convectiveHeatFlux_L
    (
        IOobject
        (
            "convectiveHeatFlux_L",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, 0, -3, 0 , 0, 0, 0),0.0)
    );

    surfaceScalarField convectiveHeatFlux_T
    (
        IOobject
        (
            "convectiveHeatFlux_T",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, 0, -3, 0 , 0, 0, 0),0.0)
    );

    volScalarField mlrSolid
    (
        IOobject
    (
            "mlrSolid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimArea/dimTime,0.0)
    );
