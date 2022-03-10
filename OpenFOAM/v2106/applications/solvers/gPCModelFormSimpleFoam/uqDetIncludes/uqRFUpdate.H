/**************************** gPCE for nut/R field ****************************/

if ( REVF )
{
    volScalarField  nutDet(turbulence->nut());
    volScalarField  tkeDet(turbulence->k());

    forAll(nut, k)
    {
        nut[k] = nutDet * expGpCoeffs[k];
        tke[k] = tkeDet * expGpCoeffs[k];
    }

    //nut[0] = nutDet;
    //tke[0] = tkeDet;
}

if ( RRSTF )
{
    volSymmTensorField  RDet(turbulence->R());
    dimensionedScalar   dimFactor("A",dimensionSet(0,2,-2,0,0,0,0), 1.0);
    volSymmTensorField  RDetDimLess(RDet/dimFactor);
    volTensorField      L(GCoeffs[0] * 0.0);

    volScalarField      A11(mag(RDetDimLess.component(0)));
    volScalarField      A21(RDetDimLess.component(1));
    volScalarField      A31(RDetDimLess.component(2));
    volScalarField      A22(mag(RDetDimLess.component(3)));
    volScalarField      A32(RDetDimLess.component(4));
    volScalarField      A33(mag(RDetDimLess.component(5)));

    // Cholesky decomposition of RDet
    L.replace(0, sqrt(A11));
    L.replace(1, A21/L.component(0));
    L.replace(2, A31/L.component(0));
    L.replace(4, sqrt(mag(A22-sqr(L.component(1)))));
    L.replace(5, (A32-L.component(2)*L.component(1))/L.component(4));
    L.replace(8, sqrt(mag(A33-sqr(L.component(2))-sqr(L.component(5)))));

    forAll(R, k)
    {
        R[k] = symm(L.T() & GCoeffs[k] & L) * dimFactor;
    }

    //R[0] = RDet;
}

