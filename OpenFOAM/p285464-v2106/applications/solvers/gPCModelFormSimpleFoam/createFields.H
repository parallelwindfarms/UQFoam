// UQ pre-processing
Info<< "--Reading UQ Inner Products calculated by chaospy\n" << endl;
#include "uqCreatefields.H"

/*********************************************************************/

// Reading fields for Deterministic solution
Info<< "\n--Creating field p\n" << endl;
volScalarField pDet
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "\n--Creating field U\n" << endl;
volVectorField UDet
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

Info<< "\n--Creating field phi\n" << endl;
surfaceScalarField phiDet
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(UDet)
);

/*********************************************************************/

word simulationType(runTime.controlDict().lookup("simulationType"));

Info<<"\n-Models for deterministic solution: " << endl;
singlePhaseTransportModel laminarTransport(UDet, phiDet);
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(UDet, phiDet, laminarTransport)
);

/*
Info<<"\n-Models for mean solution: " << endl;
singlePhaseTransportModel laminarTransport(U[0], phi[0]);
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U[0], phi[0], laminarTransport)
);
//*/
/*********************************************************************/

// Reading fields for UQ solution
PtrList<volScalarField> p(P+1);

Info<< "--Reading field p[k]\n" << endl;
forAll(p, k)
{
    Info << "Reading field p" << k << endl;
    p.set(k,new volScalarField
    (
        IOobject
        (
            "p" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ));
}

PtrList<volVectorField> U(P+1);
Info<< "\n--Reading field U[k]\n" << endl;
forAll(U, k)
{
    Info << "Reading field U" << k << endl;
    U.set(k,new volVectorField
    (
        IOobject
        (
            "U" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ));
}

PtrList<surfaceScalarField> phi(P+1);
Info<< "\n--Reading field phi[k]\n" << endl;
forAll(phi, k)
{
    Info << "Reading field phi" << k << endl;
    phi.set(k,new surfaceScalarField
    (
        IOobject
        (
            "phi" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U[k])
    ));
}

PtrList<volScalarField> nu(P+1);
Info<< "\n--Reading field nu[k]\n" << endl;
forAll(nu, k)
{
    // Change this condition if nu is a random variable/field
    if ( k == 0 )
    {
        Info << "Reading field nu" << k << endl;
        nu.set(k,new volScalarField
        (
            IOobject
            (
                "nu" + std::to_string(k),
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
           mesh
        ));
    }
    else
    {

        Info << "Reading field nu" << k << endl;
        nu.set(k,new volScalarField
        (
            IOobject
            (
                "nu" + std::to_string(k),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
           mesh,
           dimensionedScalar(sqr(dimLength)/dimTime, Zero)
        ));
    }
}

PtrList<volScalarField> nut(P+1);   // For REVF
Info<< "\n--Reading field nut[k]\n" << endl;
forAll(nut, k)
{
    // Change this condition if all the modes of nut are required in the time dirs
    Info << "Reading field nut" << k << endl;
    nut.set(k,new volScalarField
    (
        IOobject
        (
            "nut" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       turbulence->nut() * expGpCoeffs[k]
    ));
}

PtrList<volSymmTensorField> R(P+1);
Info<< "\n--Creating field R[k]\n" << endl;
forAll(R, k)
{
    // Change this condition if all the modes of nut are required in the time dirs
    Info << "Reading field R" << k << endl;
    R.set(k,new volSymmTensorField
    (
        IOobject
        (
            "R" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE  //AUTO_WRITE
        ),
        symm(U[0]*U[0])*0.0
    ));
}

PtrList<volScalarField> tke(P+1);   // FIX IT!
Info<< "\n--Reading field tke[k]\n" << endl;
forAll(tke, k)
{
    // Change this condition if all the modes of k are required in the time dirs
    if ( k == 0 )
    {
        Info << "Reading field tke" << k << endl;
        tke.set(k,new volScalarField
        (
            IOobject
            (
                "k" + std::to_string(k),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            turbulence->k() * expGpCoeffs[k]
        ));
    }
    else
    {
        Info << "Reading field tke" << k << endl;
        tke.set(k,new volScalarField
        (
            IOobject
            (
                "k" + std::to_string(k),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE  //AUTO_WRITE
            ),
            turbulence->k() * expGpCoeffs[k]
        ));
    }
}

Info<< "\n--Creating nutSigma field kSigma\n" << endl;
volScalarField tkeSigma
(
    IOobject
    (
        "kSigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    tke[1]*0.0
);


label  pRefCell = 0;
scalar pRefValue = 0.0;

//setRefCell(pDet, pimple.dict(), pRefCell, pRefValue);
setRefCell(pDet, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(pDet.name());

forAll(p, k)
{
    //setRefCell(p[k], pimple.dict(), pRefCell, pRefValue);
    setRefCell(p[k], simple.dict(), pRefCell, pRefValue);
    mesh.setFluxRequired(p[k].name());
}

#include "createFvOptions.H"

/*********************************************************************/

// UQ mean, st.dev. and their ModEff initialization (for postProcessing)


/**** Pressure Mean, Sigma ****/
/*PtrList<volScalarField> pMean(P+1);
Info<< "\n--Reading field p[k]Mean\n" << endl;
forAll(pMean, k)
{
    pMean.set(k,new volScalarField
    (
        IOobject
        (
            "pMean" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        p[k]*0.0
    ));
}*/

Info<< "\n--Creating field pSigma\n" << endl;
volScalarField pSigma
(
    IOobject
    (
        "pSigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    p[1]*0.0
);

//Info<< "\n--Creating field pMeanSigma\n" << endl;
/*volScalarField pMeanSigma
(
    IOobject
    (
        "pMeanSigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    p[1]*0.0
);*/



/**** Velocity Mean, Sigma ****/
PtrList<volVectorField> UMean(P+1);
Info<< "\n--Reading field U[k]Mean\n" << endl;
forAll(UMean, k)
{
    UMean.set(k,new volVectorField
    (
        IOobject
        (
            "UMean" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE  //AUTO_WRITE
        ),
        U[k]*0.0
    ));
}

Info<< "\n--Creating field USigma\n" << endl;
volVectorField USigma
(
    IOobject
    (
        "USigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U[1]*0.0
);

/*Info<< "\n--Creating field UMeanSigma\n" << endl;
volVectorField UMeanSigma
(
    IOobject
    (
        "UMeanSigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U[1]*0.0
);*/



/**** Turbulent viscosity Mean, Sigma ****/
/*PtrList<volScalarField> nutMean(P+1);
Info<< "\n--Reading field nut[k]Mean\n" << endl;
forAll(nutMean, k)
{
    nutMean.set(k,new volScalarField
    (
        IOobject
        (
            "nutMean" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        nut[k]*0.0
    ));
}*/

Info<< "\n--Creating nutSigma field nutSigma\n" << endl;
volScalarField nutSigma
(
    IOobject
    (
        "nutSigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    nut[1]*0.0
);

/*Info<< "\n--Creating field nutMeanSigma\n" << endl;
volScalarField nutMeanSigma
(
    IOobject
    (
        "nutMeanSigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    nut[1]*0.0
);*/



/**** Reynolds Stress Mean, Sigma ****/ // For RRSTF

/*PtrList<volSymmTensorField> RMean(P+1);
Info<< "\n--Creating field RMean[k]\n" << endl;
forAll(RMean, k)
{
    // Change this condition if all the modes of RMean are required in the time dirs
    if ( k == 0 )
    {
        RMean.set(k,new volSymmTensorField
        (
            IOobject
            (
                "RMean" + std::to_string(k),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE    //NO_WRITE
            ),
            symm(U[0]*U[0])*0.0
        ));
    }
    else
    {
        RMean.set(k,new volSymmTensorField
        (
            IOobject
            (
                "RMean" + std::to_string(k),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE    //NO_WRITE
            ),
            symm(U[0]*U[0])*0.0
        ));
    }
}*/

Info<< "\n--Creating field RSigma\n" << endl;
volSymmTensorField RSigma
(
    IOobject
    (
        "RSigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    symm(U[1]*U[1])*0.0
);

/*Info<< "\n--Creating field RMeanSigma\n" << endl;
volSymmTensorField RMeanSigma
(
    IOobject
    (
        "RMeanSigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    symm(U[0]*U[0])*0.0
);*/



/**** wallShearStress Mean, Sigma ****/
PtrList<volVectorField> tauW(P+1);
Info<< "\n--Reading field tauW[k]\n" << endl;
forAll(tauW, k)
{
    // Change this condition if all the modes of tauW are required in the time dirs
    if ( k == 0 )
    {
        tauW.set(k,new volVectorField
        (
            IOobject
            (
                "tauW" + std::to_string(k),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector(sqr(dimLength)/sqr(dimTime), Zero)
        ));
    }
    else
    {
        tauW.set(k,new volVectorField
        (
            IOobject
            (
                "tauW" + std::to_string(k),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE  //AUTO_WRITE
            ),
            mesh,
            dimensionedVector(sqr(dimLength)/sqr(dimTime), Zero)
        ));
    }
}

/*PtrList<volVectorField> tauWMean(P+1);
Info<< "\n--Reading field tauWMean[k]\n" << endl;
forAll(tauWMean, k)
{
    tauWMean.set(k,new volVectorField
    (
        IOobject
        (
            "tauWMean" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(sqr(dimLength)/sqr(dimTime), Zero)
    ));
}*/

Info<< "\n--Creating field tauWSigma\n" << endl;
volVectorField tauWSigma
(
    IOobject
    (
        "tauWSigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(sqr(dimLength)/sqr(dimTime), Zero)
);

/*Info<< "\n--Creating field tauWMeanSigma\n" << endl;
volVectorField tauWMeanSigma
(
    IOobject
    (
        "tauWMeanSigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(sqr(dimLength)/sqr(dimTime), Zero)
);*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
