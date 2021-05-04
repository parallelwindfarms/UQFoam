// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + 
    fvm::div(phi, U)
  + MRF.DDt(U)
  + fvc::div(R)
  //+ fvc::laplacian(turbulence->nut(), U, "laplacian(nuEff,U)")
  //- fvm::laplacian(turbulence->nuEff(), U)
  - fvm::laplacian(turbulence->nu(), U)
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}

/************************************************************/
