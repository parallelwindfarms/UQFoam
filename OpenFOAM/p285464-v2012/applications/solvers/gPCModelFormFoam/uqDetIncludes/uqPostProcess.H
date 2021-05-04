/********************** p, U, nut, R - Sigma **********************/
/********************** tau_w - Mean, Sigma ***********************/
/*
if (runTime.writeTime() )
{

    volScalarField     tmpNutVar (sqr(nut[0]) * 0.0);
    volScalarField     tmpPVar   (sqr(p[0])   * 0.0);
    volVectorField     tmpUVar   (cmptMultiply(U[0],U[0]) * 0.0);
    volSymmTensorField tmpRVar   (cmptMultiply(R[0],R[0]) * 0.0);
    symmTensor         Ei        (1,0,0,1,0,1);

    forAll(U, k)
    {
        R[k] *= 0.0;

        forAll(R, i)
        {
            forAll(R, j)
            {
                Cijk = M3[k][i][j];

                if ( Cijk != 0.0 )
                {
                    R[k] += - Cijk * nut[i]*dev(twoSymm(fvc::grad(U[j])));
                }
            }
        }
        R[k] += M1[k] * (2.0/3.0)*tke[k]*Ei;
        R[k] = R[k] / M1[k];

        if ( k > 0 )
        {
            tmpNutVar += M1[k] * sqr(nut[k]);
            tmpPVar   += M1[k] * sqr(p[k]);
            tmpUVar   += M1[k] * cmptMultiply(U[k], U[k]);
            tmpRVar   += M1[k] * cmptMultiply(R[k], R[k]);
        }

    }

    nutSigma = sqrt(tmpNutVar);
    pSigma   = sqrt(tmpPVar);
    for (label dir=0; dir<3; dir++)
        USigma.replace(dir, sqrt(tmpUVar.component(dir)));
    for (label dir=0; dir<6; dir++)
        RSigma.replace(dir, sqrt(tmpRVar.component(dir)));



    volSymmTensorField  Reff        (symm(U[1]*U[1]) * 0.0);
    volVectorField      tmpTauWVar  (cmptMultiply(tauW[0], tauW[0]) * 0.0);
    word                tauWPatch   (runTime.controlDict().lookup("tauWPatch"));
    label               patchID     (mesh.boundaryMesh().findPatchID(tauWPatch));

    const vectorField sf_magSf = - mesh.Sf().boundaryField()[patchID]
                               / mag(mesh.Sf().boundaryField()[patchID]);

    forAll(tauW, k)
    {
        tauW[k].boundaryFieldRef()[patchID] *= 0.0;

        forAll(tauW, i)
        {
            forAll(tauW, j)
            {
                Cijk = M3[k][i][j];

                if ( Cijk != 0.0 )
                {
                    Reff = -(nu[i]+nut[i])*dev(twoSymm(fvc::grad(U[j])));

                    tauW[k].boundaryFieldRef()[patchID] += Cijk *
                            (sf_magSf & (Reff.boundaryField()[patchID]));
                }
            }
        }

        if ( k > 0 )
        {
            tmpTauWVar += M1[k] * cmptMultiply(tauW[k], tauW[k]);
        }
    }

    for (label dir=0; dir<3; dir++)
        tauWSigma.replace(dir, sqrt(tmpTauWVar.component(dir)));

}
//*/
// End of writeTime


/****************************** Restart settings  **************************************/

label timeIndex(runTime.timeIndex());

Info << "Time index = " << timeIndex << endl;

if ( timeIndex % avgInterval == 0 )
{
    scalar dt   = runTime.deltaT().value();
    scalar time = runTime.value();
    scalar beta = 1.0/(timeIndex/avgInterval + 1.0);
    label  restartAvgTimeIndex = int(restartAvgTime/dt);

    if (timeIndex-restartAvgTimeIndex >= 0)
    {
       beta = 1.0/((timeIndex-restartAvgTimeIndex)/avgInterval + 1.0);
    }

    Info << "Time averaging at t = " << time << " beta = " << beta << endl;

    /********************** nut, p, U - Mean, Sigma, MeanSigma **********************/

//    if ( RRSTF )
//    {
//        forAll(nut, k)
//        {
//            /************ NEEDS FIXING ************/
//            nut[k] = ( (2/3 * turbulence->k() * symm(GCoeffs[k]) - R[k]) &&
//                       twoSymm(fvc::grad(U[k])) ) /
//                     ( twoSymm(fvc::grad(U[k])) && twoSymm(fvc::grad(U[k])) );
//        }
//    }

    forAll(UMean, k)
    {
//        nutMean[k]  = (1-beta)*nutMean[k] + beta*nut[k];
//        pMean[k]    = (1-beta)*pMean[k]   + beta*p[k];
        UMean[k]    = (1-beta)*UMean[k]   + beta*U[k];
    }

    volScalarField      tmpnutVar     (0.0*sqr(nut[0]));
    volScalarField      tmppVar       (0.0*sqr(p[0]));
    volVectorField      tmpUVar       (0.0*cmptMultiply(U[0], U[0]));
//    volScalarField      tmpnutMeanVar (0.0*sqr(nutMean[0]));
//    volScalarField      tmppMeanVar   (0.0*sqr(pMean[0]));
//    volVectorField      tmpUMeanVar   (0.0*cmptMultiply(UMean[0], UMean[0]));

    forAll(U, k)
    {
        if ( k > 0 )
        {
            tmpnutVar       += M1[k] * sqr(nut[k]);
            tmppVar         += M1[k] * sqr(p[k]);
            tmpUVar         += M1[k] * cmptMultiply(U[k], U[k]);
//            tmpnutMeanVar   += M1[k] * sqr(nutMean[k]);
//            tmppMeanVar     += M1[k] * sqr(pMean[k]);
//            tmpUMeanVar     += M1[k] * cmptMultiply(UMean[k], UMean[k]);
        }
    }

    nutSigma     = sqrt(tmpnutVar);
//    nutMeanSigma = sqrt(tmpnutMeanVar);

    pSigma       = sqrt(tmppVar);
//    pMeanSigma   = sqrt(tmppMeanVar);

    for (label dir=0; dir<3; dir++)
    {
        USigma.replace(dir, sqrt(tmpUVar.component(dir)));
//        UMeanSigma.replace(dir, sqrt(tmpUMeanVar.component(dir)));
    }

    /*************************** R - Mean, Sigma, MeanSigma *************************/

    if ( time > dt )
    {
        volSymmTensorField  tmpRk       (0.0*R[0]);
        volSymmTensorField  tmpRVar     (0.0*cmptMultiply(R[0], R[0]));
        //volSymmTensorField  tmpRMeanVar (0.0*cmptMultiply(RMean[0], RMean[0]));
        symmTensor          Ei          (1,0,0,1,0,1);

        if ( REVF )
        {
            forAll(R, k)
            {
                tmpRk *= 0.0;

                forAll(R, i)
                {
                    forAll(R, j)
                    {
                        Cijk = M3[k][i][j];

                        if ( Cijk != 0.0 )
                        {
                            if ( simulationType=="RAS" )
                                tmpRk += - Cijk * nut[i]*dev(twoSymm(fvc::grad(U[j])));
                            else if ( simulationType=="LES" )
                                tmpRk +=  Cijk * symm((U[i]-UMean[i])*(U[j]-UMean[j]));
                        }
                    }
                }
                if ( simulationType=="RAS" )
                    tmpRk += M1[k] * (2.0/3.0)*tke[k]*Ei;

                //RMean[k] = (1-beta)*RMean[k] + beta*tmpRk/M1[k];
                if ( k > 0)
                {
                    tmpRVar     += M1[k] * cmptMultiply(tmpRk, tmpRk);
                    //tmpRMeanVar += M1[k] * cmptMultiply(RMean[k], RMean[k]);
                }
            }
        }

        if ( RRSTF )
        {
            forAll(R, k)
            {
                //RMean[k] = (1-beta)*RMean[k] + beta*R[k];

                if ( k > 0)
                {
                    tmpRVar     += M1[k] * cmptMultiply(R[k], R[k]);
                    //tmpRMeanVar += M1[k] * cmptMultiply(RMean[k], RMean[k]);
                }
            }
        }


        for (label dir=0; dir<6; dir++)
        {
            RSigma.replace(dir, sqrt(tmpRVar.component(dir)));
            //RMeanSigma.replace(dir, sqrt(tmpRMeanVar.component(dir)));
        }
    }

    /************************** tauW - Mean, Sigma, MeanSigma ***********************/

    volSymmTensorField  Reff            (symm(U[1]*U[1]) * 0.0);
    volVectorField      tmpTauWVar      (cmptMultiply(tauW[0], tauW[0]) * 0.0);
    volVectorField      tmpTauWMeanVar  (cmptMultiply(tauW[0], tauW[0]) * 0.0);
    word                tauWPatch       (runTime.controlDict().lookup("tauWPatch"));
    label               patchID         (mesh.boundaryMesh().findPatchID(tauWPatch));

    const vectorField sf_magSf = - mesh.Sf().boundaryField()[patchID]
                               / mag(mesh.Sf().boundaryField()[patchID]);

    forAll(tauW, k)
    {
        tauW[k].boundaryFieldRef()[patchID] *= 0.0;

        forAll(tauW, i)
        {
            forAll(tauW, j)
            {
                Cijk = M3[k][i][j]/M1[k];

                if ( Cijk != 0.0 )
                {
                    Reff = -(nu[i]+nut[i])*dev(twoSymm(fvc::grad(U[j])));

                    tauW[k].boundaryFieldRef()[patchID] += Cijk *
                            (sf_magSf & (Reff.boundaryField()[patchID]));
                }
            }
        }

//        tauWMean[k] = (1-beta)*tauWMean[k] + beta*tauW[k];

        if ( k > 0 )
        {
            tmpTauWVar     += M1[k] * cmptMultiply(tauW[k], tauW[k]);
//            tmpTauWMeanVar += M1[k] * cmptMultiply(tauWMean[k], tauWMean[k]);
        }
    }

    for (label dir=0; dir<3; dir++)
    {
        tauWSigma.replace(dir, sqrt(tmpTauWVar.component(dir)));
//        tauWMeanSigma.replace(dir, sqrt(tmpTauWMeanVar.component(dir)));
    }
}
