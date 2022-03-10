#include "detCourantNo.H"

// --- Pressure-velocity PIMPLE corrector loop
while (pimple.loop())
{
    #include "detUEqn.H"

    // --- Pressure corrector loop
    while (pimple.correct())
    {
        #include "detPEqn.H"
    }

    if (pimple.turbCorr())
    {
        laminarTransport.correct();
        turbulence->correct();
    }
}
