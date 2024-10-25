/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "RVFGasSootEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "zeroGradientFvPatchFields.H"
#include "basicMultiComponentMixture.H"

#include "surfaceFields.H"

/*#include "psiCombustionModel.H"*/
/*#include "psiThermoCombustion.H"*/
#include "psiReactionThermo.H"
#include "thermoPhysicsTypes.H"
#include "eddyDissipationRVFModel.H"
#include "scalarIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(RVFGasSootEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            RVFGasSootEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::RVFGasSootEmission::RVFGasSootEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_((dict.subDict(typeName + "Coeffs"))),
    EhrrCoeff_(readScalar(coeffsDict_.lookup("EhrrCoeff"))),
    radScaling(coeffsDict_.lookupOrDefault<Switch>("radScaling",false)),
    patch2isCUP_(coeffsDict_.lookupOrDefault<Switch>("patch2isCUP",false)),
    Ehrr1_(coeffsDict_.lookupOrDefault<scalar>("Ehrr1",0.3)),
    Ehrr2_
    (
        patch2isCUP_?
        coeffsDict_.lookupOrDefault<scalar>("Ehrr2",0.368) :
        coeffsDict_.lookupOrDefault<scalar>("Ehrr2",0.3)
    ),
    patchName1_(coeffsDict_.lookup("patch1")),
    patchName2_(coeffsDict_.lookup("patch2")),
    Ehrr2CC_(coeffsDict_.lookupOrDefault<scalar>("Ehrr2CC",0.22)),
    Ehrr2PS_(coeffsDict_.lookupOrDefault<scalar>("Ehrr2PS",0.46)),
    EhrrMin_
    (
        patch2isCUP_?
        min(Ehrr1_, min(Ehrr2_, min(Ehrr2CC_, Ehrr2PS_))) :
        min(Ehrr1_, Ehrr2_)
    ),
    EhrrMax_
    (
        patch2isCUP_?
        max(Ehrr1_, max(Ehrr2_, max(Ehrr2CC_, Ehrr2PS_))) :
        max(Ehrr1_, Ehrr2_)
    ),
    tstartXrLocal_(coeffsDict_.lookupOrDefault<scalar>("localRadFracStartTime",40.))
{
    if (patch2isCUP_)
    {
        Info << "patch2 radiant fraction will be mlr-weighted based on CUP pyrolysis model"<<nl 
             << "    Ehrr2CC = " << Ehrr2CC_ <<nl
             << "    Ehrr2PS = " << Ehrr2PS_ <<nl
             << nl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::RVFGasSootEmission::~RVFGasSootEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::RVFGasSootEmission::aCont(const label bandI) const
{

    tmp<volScalarField> a
    (
        new volScalarField
        (
            IOobject
            (
                "aCont" + name(bandI),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0),
            zeroGradientFvPatchVectorField::typeName
        )
    );

    return a;

}


Foam::tmp<Foam::volScalarField>
Foam::radiation::RVFGasSootEmission::eCont(const label bandI) const
{
    tmp<volScalarField> e
    (
        new volScalarField
        (
            IOobject
            (
                "eCont" + name(bandI),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("e", dimless/dimLength, 0.0)
        )
    );

    return e;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::RVFGasSootEmission::ECont(const label bandI) const
{
    tmp<volScalarField> E
    (
        new volScalarField
        (
            IOobject
            (
                "ECont" + name(bandI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );

    scalar RadFraction = 0;

    if (radScaling)
    {

        //TODO: this doesn't need to be recomputed for each ILambda solve

        const surfaceScalarField& phi = mesh_.lookupObject<surfaceScalarField>("phi");

        scalar mlr1(0.0);

        // kvm added this ....

        forAll(patchName1_,i)
        {
            const label patchI = mesh_.boundaryMesh().findPatchID(patchName1_[i]);
            if(patchI<0)
            {
                FatalErrorIn("radScaling.H")
                    << "patch " << patchName1_[i] << " not found" << nl
                    << abort(FatalError);
            }
            mlr1 += -gSum(phi.boundaryField()[patchI]);
        }

        scalar mlr2(0.0);
        scalar mlr2xEhrr2(0.0);

        forAll(patchName2_,i)
        {
            const label patchI = mesh_.boundaryMesh().findPatchID(patchName2_[i]);
            if(patchI<0)
            {
                FatalErrorIn("radScaling.H")
                    << "patch " << patchName2_[i] << " not found" << nl
                    << abort(FatalError);
            }

            if (patch2isCUP_)
            {
                const scalarIOList& mlrPS = mesh_.lookupObject<scalarIOList>("mlrPS");
                const scalarIOList& mlrCC = mesh_.lookupObject<scalarIOList>("mlrCC");
                const scalar& mlr2CC = mlrCC[patchI];
                const scalar& mlr2PS = mlrPS[patchI];

                mlr2 += (mlr2CC + mlr2PS);
                mlr2xEhrr2 += (Ehrr2CC_*mlr2CC + Ehrr2PS_*mlr2PS);
            }
            else
            {
                scalar mlr2p = -gSum(phi.boundaryField()[patchI]);
                mlr2 += mlr2p;
                mlr2xEhrr2 += mlr2p*Ehrr2_;
            }

        }

        if(debug)
        {
            Info << "mlr for patches " << patchName1_ << " is " << mlr1 << endl;
            Info << "mlr for patches " << patchName2_ << " is " << mlr2 << endl;
        }

        RadFraction = (mlr1*Ehrr1_ + mlr2xEhrr2)
                    / max(SMALL, (mlr1 + mlr2));
        RadFraction = max(EhrrMin_,RadFraction);
        //debug Info << "RadFraction " << RadFraction << endl;
    }
    else
    {
        RadFraction = EhrrCoeff_;
    }


    if (mesh_.foundObject<volScalarField>("Qdot"))
    {
        const volScalarField& Qdot =
            mesh_.lookupObject<volScalarField>("Qdot");
        if (Qdot.dimensions() == dimEnergy/dimTime/dimVolume)
        {
            //const bool& RVFModelActivated = mesh_.lookupObject<bool>("RVFModelActivated_");
            //const bool& fixedXr = mesh_.lookupObject<bool>("fixedXr_");
            //const bool fixedXr = this->mesh().template lookupObject<bool>("fixedXr_");
            //const combustionModels::psiCombustionModel& combM = mesh_.lookupObject<combustionModels::psiCombustionModel>("eddyDissipationRVFModel");
            //if(combM.modelActivated() && !combM.fixedXr())
            //if(mesh().time().value()>40)

            /*typedef combustionModels::psiCombustionModel modelType;*/
            typedef combustionModels::singleStepCombustion<psiReactionThermo,gasHThermoPhysics> modelType;
            //const modelType& combustionModel = mesh_.lookupObject<modelType>("eddyDissipationRVFModel");
            const modelType& combustionModel = mesh_.lookupObject<modelType>("combustionProperties");
            const combustionModels::eddyDissipationRVFModel<psiReactionThermo,gasHThermoPhysics>& combM = dynamic_cast<const combustionModels::eddyDissipationRVFModel<psiReactionThermo,gasHThermoPhysics>&>(combustionModel);
            const scalar& currentT(mesh().time().value());
            const Switch variableXr = 
                         ( 
		             patch2isCUP_? 
			     ( combM.variableXr() && (currentT > tstartXrLocal_) ) :
			     combM.variableXr()
			 );
            if (variableXr)
            {
                //- Ning: using local radiant fraction from simple soot model
                const volScalarField& XrLocal =
                    mesh_.lookupObject<volScalarField>("XrFlame");
                //E.ref().ref() = XrLocal*Qdot;
                const scalar zero = 0.0;
                const scalar one = 1.0;
                if (patch2isCUP_)
                {
                    scalar blendR(max(zero,(min(one,(currentT-tstartXrLocal_)/10.0))));
                    const volScalarField extFactor(blendR*max(zero, min(one, (one - XrLocal/Ehrr2_))));
                    E.ref().ref() = (1.0 - extFactor)*RadFraction*Qdot;
                }
                else
                {
                    scalar blendR(max(zero,(min(one,(currentT-combM.RVFactivationTime())/10.0))));
                    E.ref().ref() = (blendR*XrLocal+(1.0-blendR)*RadFraction)*Qdot;
                }
            }
            else
            {
                E.ref().ref() = RadFraction*Qdot;

                static word timeName = "null";
                if (timeName != mesh().time().timeName())
                {
                    Info << "Radiant Fraction is " << RadFraction << endl;
                    timeName = mesh().time().timeName();
                }
            }
        }
        else
        {
            Info << "Qdot dimensions incorrect" << endl;
        }
    }

    return E;
}


// ************************************************************************* //
