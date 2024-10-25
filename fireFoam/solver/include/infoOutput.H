 //- Convective heat flux
 convectiveHeatFlux_T = 
     - fvc::interpolate(turbulence->alphaEff())*fvc::interpolate(thermo.Cp())*fvc::snGrad(T);        

     convectiveHeatFlux_L = 
         - fvc::interpolate(turbulence->alpha())*fvc::interpolate(thermo.Cp())*fvc::snGrad(T);


//- Solid mass loss rate
forAll(mlrSolid.boundaryField(), patchi)
{
    if (mesh.boundary()[patchi].type() == "mappedWall")
    {
            mlrSolid.boundaryFieldRef()[patchi] =
            (
                - phi.boundaryField()[patchi]
                /mesh.boundary()[patchi].magSf()
            );
    }
}

//- For flame visualization
if (thermo.composition().contains("O2")) 
{
    ft = (fu*s-O2+YO2Inf)/(s+YO2Inf);
}
else
{
    ft = fu;
}

ft.max(0.0);
ft.min(1.0);

