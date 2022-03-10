// --- Explicit iteration cycles b/w systems of Uk-pk
#include "U0CourantNo.H"

while (expItr < expItrMax)
{
    expItr++;

    Info<< nl << "expItr: cycle " << expItr << endl;

    // --- Loop over velocity modes (reverse)
    forAll(U, k)
    {
        // --- Efficiently updating UQ modes
        if(k<=Ptrunc)
        {
            Info<< "gPC: mode " << k << endl;

            // --- Pressure-velocity PIMPLE/SIMPLE corrector loop

              #include "UkCourantNo.H"

              // --- Momentum predictor
              #include "uqUEqn.H"

              // --- Pressure corrector loop
              #include "uqPEqn.H"

            Info<< "------------------------------------" << endl;

        }

    } // --- End of velocity modes loop

} // --- End of explicit cycles
