// Solve the Momentum equation

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(UDet) + fvm::div(phiDet, UDet)
  + turbulence->divDevReff(UDet)
 ==
    fvOptions(UDet)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(pDet));

    fvOptions.correct(UDet);
}
