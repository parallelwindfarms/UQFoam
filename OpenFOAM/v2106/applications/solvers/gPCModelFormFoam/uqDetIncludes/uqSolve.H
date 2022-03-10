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
        //if((expItrMax-expItr >= (P-k) || expItr <= k+1) && (k<=Ptrunc))
        if(k<=Ptrunc)
        {
            Info<< "gPC: mode " << k << endl;

            // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
                  #include "UkCourantNo.H"

                  // --- Momentum predictor
                  #include "uqUEqn.H"

                  // --- Pressure corrector loop
                  while (pimple.correct())
                  {
                      #include "uqPEqn.H"
                  }

            } // --- End of PIMPLE loop

            Info<< "------------------------------------" << endl;

        }

    } // --- End of velocity modes loop

} // --- End of explicit cycles
