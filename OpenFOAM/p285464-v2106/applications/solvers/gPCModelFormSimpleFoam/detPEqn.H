volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), UDet, pDet));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
);

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(pDet)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(pDet);
}

if (pimple.nCorrPISO() <= 1)
{
    //tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(pDet, UDet, phiHbyA, rAtU());

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), pDet) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(pDet.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phiDet = phiHbyA - pEqn.flux();
    }
}

//#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
pDet.relax();

UDet = HbyA - rAtU*fvc::grad(pDet);
UDet.correctBoundaryConditions();
fvOptions.correct(UDet);

// Correct Uf if the mesh is moving
//fvc::correctUf(UfDet, UDet, phiDet);

// Make the fluxes relative to the mesh motion
//fvc::makeRelative(phiDet, UDet);
