
IOdictionary uqInfo
(
    IOobject
    (
        "uqInfo",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

/***************************************************************/

Info<< "--Reading type of model form UQ" << endl;
word randField(uqInfo.lookup("randField"));
bool REVF {false}, RRSTF {false};

if ( randField=="REVF" )
{
    REVF = true;
    Info<< "Random eddy viscosity field (REV) used for UQ" << endl;
}

if ( randField=="RRSTF" )
{
    RRSTF = true;
    Info<< "Random Reynolds stress tensor (RRSTF) field used for UQ" << endl;
}

Info<< "--Reading total number of modes (P)" << endl;
label P(readLabel(uqInfo.lookup("P")));
Info<< "\n  P = " << P << " => k = [0," << P << "]\n" << endl;

label  Pplus1 = P+1;

Info<< "--Reading totla number of dimensions (d)" << endl;
label d(readLabel(uqInfo.lookup("d")));
Info<< "\n  d = " << d << endl;

/***************************************************************/

PtrList<volScalarField> expGpCoeffs(P+1);
Info<< "\n--Reading nut gPCE coefficients expGpCoeffs[k]\n" << endl;
forAll(expGpCoeffs, k)
{
    Info << "Reading field nutCoeffs" << k << endl;
    expGpCoeffs.set(k,new volScalarField
    (
        IOobject
        (
            "expGpCoeffs" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            ((REVF)? (IOobject::AUTO_WRITE):(IOobject::NO_WRITE))
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    ));
}

PtrList<volTensorField> GCoeffs(P+1);
Info<< "\n--Reading G gPCE coefficients GCoeffs[k]\n" << endl;
forAll(GCoeffs, k)
{
    Info << "Reading field GCoeffs" << k << endl;
    GCoeffs.set(k,new volTensorField
    (
        IOobject
        (
            "GCoeffs" + std::to_string(k),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            ((RRSTF)? (IOobject::AUTO_WRITE):(IOobject::NO_WRITE))
        ),
        mesh,
        dimensionedTensor(dimless, Zero)
    ));
}

/***************************************************************/

scalar Ckk(0.0);
scalar M1[Pplus1];

for (label k=0; k< Pplus1; k++)
  M1[k] = uqInfo.lookupOrDefault<scalar>
               ("M["+std::to_string(k)+"]",0.0);


scalar Cijk(0.0);
scalar M3[Pplus1][Pplus1][Pplus1];

for (label i=0; i < Pplus1; i++)
  for (label j=0; j< Pplus1; j++)
    for (label k=0; k< Pplus1; k++)
      M3[i][j][k] = uqInfo.lookupOrDefault<scalar>
                   ("M["+std::to_string(i)+
                    "]["+std::to_string(j)+
                    "]["+std::to_string(k)+"]",0.0);

