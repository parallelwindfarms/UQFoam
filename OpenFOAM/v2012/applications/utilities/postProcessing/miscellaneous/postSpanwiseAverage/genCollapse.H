    forAll(fieldNames, i)
    {
        //Create the IOobject for the field 
        IOobject fieldHeader
        (
            fieldNames[i],
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        //set the path
        fileName path(fieldHeader.rootPath()/fieldHeader.caseName()
        			  /"postProcessing"/"collapsedFields"/fieldHeader.instance());
    	mkDir(path);

        //Check the type of the field
        word type = "noField";
		
		if	(	
				std::regex_match (fieldNames[i], std::regex("(p)(.*)"))		||	
				std::regex_match (fieldNames[i], std::regex("(nut)(.*)"))	||
				std::regex_match (fieldNames[i], std::regex("(yPl)(.*)"))
			)
			type = "volScalarField";
			
		else if	
			(	
				std::regex_match (fieldNames[i], std::regex("(UPrime)(.*)"))	||	
				std::regex_match (fieldNames[i], std::regex("(U0Prime)(.*)"))	||
				std::regex_match (fieldNames[i], std::regex("(R)(.*)"))	
			)
			type = "volSymmTensorField";
			
		else if 
			(	
				std::regex_match (fieldNames[i], std::regex("(U)(.*)"))		||
				std::regex_match (fieldNames[i], std::regex("(tau)(.*)"))	
			)
			type = "volVectorField";
	
		// Averaging fields
		
		if ( type == "volScalarField" )
        {
            //Check if the field was found, if not, warn and do nothing
		    if (!fieldHeader.typeHeaderOk<volScalarField>(true))
		    {
		        Info<< "    No " << fieldNames[i] << " field" << endl;
		        continue;
		    }            
            
            volScalarField field
            (
                fieldHeader,
                mesh
            );
            
			volScalarField avgScalar
			(
			IOobject
			(
				fieldNames[i]+"SpAvg",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
				mesh,
				dimensionedScalar
				(
					fieldNames[i]+"SpAvg",
					dimPressure, 
					0.0
				)
			);
			
			scalarField spanavgScalar(meshIndexing.collapse(field));
    
    		avgScalar.boundaryFieldRef() = field.boundaryFieldRef();
    		
    		forAll(avgScalar.ref(), celli)
			{
			   avgScalar.ref()[celli] = spanavgScalar[celli];
			}
			
			avgScalar.write();

        }
        
        
		else if ( type == "volVectorField" )
        {
            //Check if the field was found, if not, warn and do nothing
		    if (!fieldHeader.typeHeaderOk<volVectorField>(true))
		    {
		        Info<< "    No " << fieldNames[i] << " field" << endl;
		        continue;
		    }            
            
            volVectorField field
            (
                fieldHeader,
                mesh
            );

			volScalarField avgVecX
			(
				IOobject
				(
					"avgVecX",
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
					mesh,
					dimensionedScalar("avgVecX", dimVelocity, 0.0)
			);
		 
			volScalarField avgVecY
			(
				IOobject
				(
					"avgVecY",
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
					mesh,
					dimensionedScalar("avgVecY", dimVelocity, 0.0)
			);

			volScalarField avgVecZ
			(
				IOobject
				(
					"avgVecZ",
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
					mesh,
					dimensionedScalar("avgVecZ", dimVelocity, 0.0)
			);
			
			volVectorField avgVec
			(
				IOobject
				(
					fieldNames[i]+"SpAvg",
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
					mesh,
					dimensionedVector
					(
						    fieldNames[i]+"SpAvg",
						    dimVelocity,
						    vector::zero          
					)        
			);
			
			scalarField spanavgVecX(meshIndexing.collapse(field.component(vector::X)()));
			scalarField spanavgVecY(meshIndexing.collapse(field.component(vector::Y)()));
			scalarField spanavgVecZ(meshIndexing.collapse(field.component(vector::Z)()));
			
			avgVecX.boundaryFieldRef() = field.boundaryFieldRef().component(vector::X);
			avgVecY.boundaryFieldRef() = field.boundaryFieldRef().component(vector::Y);
			avgVecZ.boundaryFieldRef() = field.boundaryFieldRef().component(vector::Z);

			forAll(avgVecX.ref(), celli)
			{
			   avgVecX.ref()[celli] = spanavgVecX[celli];
			   avgVecY.ref()[celli] = spanavgVecY[celli];
			   avgVecZ.ref()[celli] = spanavgVecZ[celli];
			}
		 
			avgVec.replace(0,avgVecX); 
			avgVec.replace(1,avgVecY);
			avgVec.replace(2,avgVecZ);

			avgVec.write();
 
        }
        
        
        
		else if ( type == "volSymmTensorField" )
        {
            //Check if the field was found, if not, warn and do nothing
		    if (!fieldHeader.typeHeaderOk<volSymmTensorField>(true))
		    {
		        Info<< "    No " << fieldNames[i] << " field" << endl;
		        continue;
		    }            
            
            volSymmTensorField field
            (
                fieldHeader,
                mesh
            );
            
			volScalarField avgTensorXX
			(
			IOobject
			(
				"avgTensorXX",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
				mesh,
				dimensionedScalar("avgTensorXX", dimVelocity*dimVelocity, 0.0)
			);
			volScalarField avgTensorYY
			(
			IOobject
			(
				"avgTensorYY",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
				mesh,
				dimensionedScalar("avgTensorYY", dimVelocity*dimVelocity, 0.0)
			);
			
			volScalarField avgTensorZZ
			(
			IOobject
			(
				"avgTensorZZ",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
				mesh,
				dimensionedScalar("avgTensorZZ", dimVelocity*dimVelocity, 0.0)
			);
			volScalarField avgTensorXY
			(
				IOobject
				(
				    "avgTensorXY",
				    runTime.timeName(),
				    mesh,
				    IOobject::NO_READ,
				    IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("avgTensorXY", dimVelocity*dimVelocity, 0.0)
			);
			volScalarField avgTensorXZ
			(
				IOobject
				(
				    "avgTensorXZ",
				    runTime.timeName(),
				    mesh,
				    IOobject::NO_READ,
				    IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("avgTensorXZ", dimVelocity*dimVelocity, 0.0)
			);
			volScalarField avgTensorYZ
			(
				IOobject
				(
				    "avgTensorYZ",
				    runTime.timeName(),
				    mesh,
				    IOobject::NO_READ,
				    IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("avgTensorYZ", dimVelocity*dimVelocity, 0.0)
			);
			
			volSymmTensorField avgTensor
			(
				IOobject
				(
					fieldNames[i]+"SpAvg",
				    runTime.timeName(),
				    mesh,
				    IOobject::NO_READ,
				    IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedSymmTensor
				(
						fieldNames[i]+"SpAvg",
				        dimVelocity*dimVelocity,
				        symmTensor::zero
				)
			);
            
            scalarField spanavgTensorXX(meshIndexing.collapse(field.component(symmTensor::XX)()));
            scalarField spanavgTensorYY(meshIndexing.collapse(field.component(symmTensor::YY)()));	
            scalarField spanavgTensorZZ(meshIndexing.collapse(field.component(symmTensor::ZZ)()));		
            scalarField spanavgTensorXY(meshIndexing.collapse(field.component(symmTensor::XY)()));
            scalarField spanavgTensorXZ(meshIndexing.collapse(field.component(symmTensor::XZ)()));            
            scalarField spanavgTensorYZ(meshIndexing.collapse(field.component(symmTensor::YZ)()));            

   			avgTensorXX.boundaryFieldRef() = field.boundaryFieldRef().component(symmTensor::XX);
			avgTensorYY.boundaryFieldRef() = field.boundaryFieldRef().component(symmTensor::YY);
			avgTensorZZ.boundaryFieldRef() = field.boundaryFieldRef().component(symmTensor::ZZ);
			avgTensorXY.boundaryFieldRef() = field.boundaryFieldRef().component(symmTensor::XY);
			avgTensorXZ.boundaryFieldRef() = field.boundaryFieldRef().component(symmTensor::XZ);
			avgTensorYZ.boundaryFieldRef() = field.boundaryFieldRef().component(symmTensor::YZ);
			
			forAll(avgTensorXX.ref(), celli)
			{
			   avgTensorXX.ref()[celli] = spanavgTensorXX[celli];
			   avgTensorYY.ref()[celli] = spanavgTensorYY[celli];
			   avgTensorZZ.ref()[celli] = spanavgTensorZZ[celli];
			   avgTensorXY.ref()[celli] = spanavgTensorXY[celli];
			   avgTensorXZ.ref()[celli] = spanavgTensorXZ[celli];
			   avgTensorYZ.ref()[celli] = spanavgTensorYZ[celli];
			}
		  
		 
			avgTensor.replace(0, avgTensorXX);
			avgTensor.replace(3, avgTensorYY);
			avgTensor.replace(5, avgTensorZZ);
			avgTensor.replace(1, avgTensorXY);
			avgTensor.replace(2, avgTensorXZ);
			avgTensor.replace(4, avgTensorYZ);

			avgTensor.write();   
        }
        
        
		else
        {
            Info << "    Unsupported field, ignoring field " << fieldNames[i] << endl;
            Info << "    Please add this field in collapse.H " << endl;
        }
	
	}



