//#include "heardersForQtParsing.H"

// --- Solve the Momentum equation
Ckk = M1[k];

fvVectorMatrix UEqn (- Ckk * fvOptions(U[k]));

if ( REVF )
{
    forAll(U, j)
    {
        forAll(U, i)
        {
            Cijk = M3[k][j][i];

            if ( Cijk != 0.0 )
            {
                if ( j==k )
                {
                    UEqn += Cijk * fvm::div(phi[i], U[j]);
                    UEqn -= Cijk * fvm::laplacian(nut[i]+nu[i], U[j]);
                }
                else
                {
                    UEqn += Cijk * fvc::div(phi[i], U[j]);
                    UEqn -= Cijk * fvc::laplacian(nut[i]+nu[i], U[j]);
                }
                //UEqn -= Cijk * fvc::div((nut[i]+nu[i])*dev2(T(fvc::grad(U[j]))));
            }
        }
    }
}

if ( RRSTF || eigPertRRSTF )
{
    forAll(U, j)
    {
        forAll(U, i)
        {
            Cijk = M3[k][j][i];

            if ( Cijk != 0.0 )
            {
                if ( j==k )
                {
                    UEqn += Cijk * fvm::div(phi[i], U[j]);
                    UEqn -= Cijk * fvm::laplacian(nu[i], U[j]);
                }
                else
                {
                    UEqn += Cijk * fvc::div(phi[i], U[j]);
                    UEqn -= Cijk * fvc::laplacian(nu[i], U[j]);
                }
                //UEqn -= Cijk * fvc::div((nu[i])*dev2(T(fvc::grad(U[j]))));
            }
        }
    }

    UEqn += Ckk * fvc::div(R[k]);
}

UEqn *= 1.0/Ckk;

UEqn.relax();

fvOptions.constrain(UEqn);

if (simple.momentumPredictor())
{
    solve(UEqn == - fvc::grad(p[k]));

    fvOptions.correct(U[k]);
}

