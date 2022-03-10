//     fileName path(UMean.rootPath()/UMean.caseName()/"span-averaged_field");
//     mkDir(path);

// VAR1



    volScalarField avgUx
    (
	IOobject
	(
	    "avgUx",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("avgUx", dimVelocity, 0.0)
    );
 
    volScalarField avgUy
    (
	IOobject
	(
	    "avgUy",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("avgUy", dimVelocity, 0.0)
    );

    volScalarField avgUz
    (
	IOobject
	(
	    "avgUz",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("avgUz", dimVelocity, 0.0)
    );

    volScalarField avgP
    (
	IOobject
	(
	    "avgP",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("avgP", dimPressure, 0.0)
    );


    volScalarField pP2M
    (
	IOobject
	(
	    "pP2M",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("pP2M", dimPressure*dimPressure, 0.0)
    );

    volScalarField UPrime2MeanXX
    (
	IOobject
	(
	    "UPrime2MeanXX",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("UPrime2MeanXX", dimVelocity*dimVelocity, 0.0)
    );
    volScalarField UPrime2MeanYY
    (
	IOobject
	(
	    "UPrime2MeanYY",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("UPrime2MeanYY", dimVelocity*dimVelocity, 0.0)
    );
    
    volScalarField UPrime2MeanZZ
    (
	IOobject
	(
	    "UPrime2MeanZZ",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedScalar("UPrime2MeanZZ", dimVelocity*dimVelocity, 0.0)
    );
    volScalarField UPrime2MeanXY
    (
        IOobject
        (
            "UPrime2MeanXY",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("UPrime2MeanXY", dimVelocity*dimVelocity, 0.0)
    );
    volScalarField UPrime2MeanXZ
    (
        IOobject
        (
            "UPrime2MeanXZ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("UPrime2MeanXZ", dimVelocity*dimVelocity, 0.0)
    );
    volScalarField UPrime2MeanYZ
    (
        IOobject
        (
            "UPrime2MeanYZ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("UPrime2MeanYZ", dimVelocity*dimVelocity, 0.0)
    );
    
    volVectorField Uvec
    (
	IOobject
	(
	    "Uvec",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
        mesh,
        dimensionedVector
        (
                "Uvec",
                dimVelocity,
                vector::zero          
        )        
    );
    
    volSymmTensorField UP2M
    (
        IOobject
        (
            "UP2M",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor
        (
                "UP2M",
                dimVelocity*dimVelocity,
                symmTensor::zero
        )
    );

//VAR1

    scalarField spanavgUMeanX
    (
        meshIndexing.collapse(UMean.component(vector::X)())
    );

    avgUx.boundaryFieldRef() = UMean.boundaryFieldRef().component(vector::X);

 //   avgUx.write();    
// VAR2
    scalarField spanavgUMeanY
    (
        meshIndexing.collapse(UMean.component(vector::Y)())
    );
    
    avgUy.boundaryFieldRef() = UMean.boundaryFieldRef().component(vector::Y);

//    avgUy.write();
// VAR3
    scalarField spanavgUMeanZ
    (
        meshIndexing.collapse(UMean.component(vector::Z)())
    );
    
    avgUz.boundaryFieldRef() = UMean.boundaryFieldRef().component(vector::Z);

//    avgUz.write();
// VAR4
    scalarField spanavgpMean
    (
        meshIndexing.collapse(pMean)
    );
    
    avgP.boundaryFieldRef() = pMean.boundaryFieldRef();

    scalarField spanavgpPrime2
    (
        meshIndexing.collapse(pPrime2Mean)
    );
    
    pP2M.boundaryFieldRef() = pPrime2Mean.boundaryFieldRef();

    scalarField spanavgUPrim2XX
    (
        meshIndexing.collapse(UPrime2Mean.component(symmTensor::XX)())
    );
    
    UPrime2MeanXX.boundaryFieldRef() = UPrime2Mean.boundaryFieldRef().component(symmTensor::XX);
//    UPrime2MeanXX.write();

    

    scalarField spanavgUPrim2YY
    (
        meshIndexing.collapse(UPrime2Mean.component(symmTensor::YY)())
    );
    
    UPrime2MeanYY.boundaryFieldRef() = UPrime2Mean.boundaryFieldRef().component(symmTensor::YY);
//    UPrime2MeanYY.write();
    
    
    
    scalarField spanavgUPrim2ZZ
    (
        meshIndexing.collapse(UPrime2Mean.component(symmTensor::ZZ)())
    );
    
    UPrime2MeanZZ.boundaryFieldRef() = UPrime2Mean.boundaryFieldRef().component(symmTensor::ZZ);
//    UPrime2MeanZZ.write();
    scalarField spanavgUPrim2XY
    (
        meshIndexing.collapse(UPrime2Mean.component(symmTensor::XY)())
    );

    UPrime2MeanXY.boundaryFieldRef() = UPrime2Mean.boundaryFieldRef().component(symmTensor::XY);

    scalarField spanavgUPrim2XZ
    (
        meshIndexing.collapse(UPrime2Mean.component(symmTensor::XZ)())
    );

    UPrime2MeanXZ.boundaryFieldRef() = UPrime2Mean.boundaryFieldRef().component(symmTensor::XZ);

    scalarField spanavgUPrim2YZ
    (
        meshIndexing.collapse(UPrime2Mean.component(symmTensor::YZ)())
    );

    UPrime2MeanYZ.boundaryFieldRef() = UPrime2Mean.boundaryFieldRef().component(symmTensor::YZ);

     forAll(avgUx.ref(), celli)
    {
       avgUx.ref()[celli] = spanavgUMeanX[celli];
       avgUy.ref()[celli] = spanavgUMeanY[celli];
       avgUz.ref()[celli] = spanavgUMeanZ[celli];
       UPrime2MeanXX.ref()[celli] = spanavgUPrim2XX[celli];
       UPrime2MeanYY.ref()[celli] = spanavgUPrim2YY[celli];
       UPrime2MeanZZ.ref()[celli] = spanavgUPrim2ZZ[celli];
       UPrime2MeanXY.ref()[celli] = spanavgUPrim2XY[celli];
       UPrime2MeanXZ.ref()[celli] = spanavgUPrim2XZ[celli];
       UPrime2MeanYZ.ref()[celli] = spanavgUPrim2YZ[celli];
       avgP.ref()[celli] = spanavgpMean[celli];
       pP2M.ref()[celli] = spanavgpPrime2[celli];

    }
  
 
   
    Uvec.replace(0,avgUx); 
    Uvec.replace(1,avgUy);
    Uvec.replace(2,avgUz);
    UP2M.replace(0, UPrime2MeanXX);
    UP2M.replace(3, UPrime2MeanYY);
    UP2M.replace(5, UPrime2MeanZZ);
    UP2M.replace(1, UPrime2MeanXY);
    UP2M.replace(2, UPrime2MeanXZ);
    UP2M.replace(4, UPrime2MeanYZ);

//pP2M.write();         
avgP.write();
Uvec.write();
UP2M.write();        
 
