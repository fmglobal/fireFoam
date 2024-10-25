Info<< "\nConstructing reacting cloud" << endl;
basicReactingCloud suppressionParcels
(
    "suppressionCloud",
    rho,
    U,
    g,
    slgThermo
);
basicReactingCloud fuelParcels
(
    "fuelCloud",
    rho,
    U,
    g,
    slgThermo
);
