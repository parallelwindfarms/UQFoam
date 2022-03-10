    //Get the corrdinates
    const scalarField& y = channelIndexing.y();
    scalarField allY(y.size()+2);

    allY[0] = yBot;
    forAll(y, yI)
    {
        allY[yI+1] = y[yI];
    }
    allY[allY.size()-1] = yTop;

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
				std::regex_match (fieldNames[i], std::regex("(UPrime)(.*)"))		||	
				std::regex_match (fieldNames[i], std::regex("(U0Prime)(.*)"))	||
				std::regex_match (fieldNames[i], std::regex("(R)(.*)"))	
			)
			type = "volSymmTensorField";
			
		else if 
			(	
				std::regex_match (fieldNames[i], std::regex("(U)(.*)"))
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

            scalarField values = channelIndexing.collapse(field);

            scalar valBot = patchAverage<volScalarField, scalar>(mesh, fieldHeader,
                                                         areaBot, patchBot);
            scalar valTop = patchAverage<volScalarField, scalar>(mesh, fieldHeader,
                                                         areaTop, patchTop);
            scalarField allValues(values.size()+2);

            allValues[0] = valBot;
            forAll(values, valI)
            {
                allValues[valI+1] = values[valI];
            }
            allValues[allValues.size()-1] = valTop;

            makeGraph(allY, allValues, fieldNames[i], path, gFormat);
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

            scalarField valuesX = channelIndexing.collapse(field.component(vector::X)());
            scalarField valuesY = channelIndexing.collapse(field.component(vector::Y)());
            scalarField valuesZ = channelIndexing.collapse(field.component(vector::Z)());

            vector valBot = patchAverage<volVectorField, vector>(mesh, fieldHeader,
                                                         areaBot, patchBot);
            vector valTop = patchAverage<volVectorField, vector>(mesh, fieldHeader,
                                                         areaTop, patchTop);
            scalarField allValuesX(valuesX.size()+2);
            scalarField allValuesY(valuesY.size()+2);
            scalarField allValuesZ(valuesZ.size()+2);

            allValuesX[0] = valBot.x();
            allValuesY[0] = valBot.y();
            allValuesZ[0] = valBot.z();

            forAll(valuesX, valI)
            {
                allValuesX[valI+1] = valuesX[valI];
                allValuesY[valI+1] = valuesY[valI];
                allValuesZ[valI+1] = valuesZ[valI];
            }
            allValuesX[allValuesX.size()-1] = valTop.x();
            allValuesY[allValuesY.size()-1] = valTop.y();
            allValuesZ[allValuesZ.size()-1] = valTop.z();

            makeGraph(allY, allValuesX, fieldNames[i]+"_X", path, gFormat);
            makeGraph(allY, allValuesY, fieldNames[i]+"_Y", path, gFormat);
            makeGraph(allY, allValuesZ, fieldNames[i]+"_Z", path, gFormat);
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

            scalarField valuesXX = channelIndexing.collapse(field.component(symmTensor::XX)());
            scalarField valuesXY = channelIndexing.collapse(field.component(symmTensor::XY)());
            scalarField valuesXZ = channelIndexing.collapse(field.component(symmTensor::XZ)());
            scalarField valuesYY = channelIndexing.collapse(field.component(symmTensor::YY)());
            scalarField valuesYZ = channelIndexing.collapse(field.component(symmTensor::YZ)());
            scalarField valuesZZ = channelIndexing.collapse(field.component(symmTensor::ZZ)());

            symmTensor valBot = patchAverage<volSymmTensorField, symmTensor>(mesh, fieldHeader,
                                                         areaBot, patchBot);
            symmTensor valTop = patchAverage<volSymmTensorField, symmTensor>(mesh, fieldHeader,
                                                         areaTop, patchTop);

            scalarField allValuesXX(valuesXX.size()+2);
            scalarField allValuesXY(valuesXY.size()+2);
            scalarField allValuesXZ(valuesXZ.size()+2);
            scalarField allValuesYY(valuesXX.size()+2);
            scalarField allValuesYZ(valuesXY.size()+2);
            scalarField allValuesZZ(valuesXZ.size()+2);

            allValuesXX[0] = valBot.component(symmTensor::XX);
            allValuesXY[0] = valBot.component(symmTensor::XY);
            allValuesXZ[0] = valBot.component(symmTensor::XZ);
            allValuesYY[0] = valBot.component(symmTensor::YY);
            allValuesYZ[0] = valBot.component(symmTensor::YZ);
            allValuesZZ[0] = valBot.component(symmTensor::ZZ);

            forAll(valuesXX, valI)
            {
                allValuesXX[valI+1] = valuesXX[valI];
                allValuesXY[valI+1] = valuesXY[valI];
                allValuesXZ[valI+1] = valuesXZ[valI];
                allValuesYY[valI+1] = valuesYY[valI];
                allValuesYZ[valI+1] = valuesYZ[valI];
                allValuesZZ[valI+1] = valuesZZ[valI];
            }

            allValuesXX[allValuesXX.size()-1] = valTop.component(symmTensor::XX);
            allValuesXY[allValuesXY.size()-1] = valTop.component(symmTensor::XY);
            allValuesXZ[allValuesXZ.size()-1] = valTop.component(symmTensor::XZ);
            allValuesYY[allValuesYY.size()-1] = valTop.component(symmTensor::YY);
            allValuesYZ[allValuesYZ.size()-1] = valTop.component(symmTensor::YZ);
            allValuesZZ[allValuesZZ.size()-1] = valTop.component(symmTensor::ZZ);

            makeGraph(allY, allValuesXX, fieldNames[i]+"_XX", path, gFormat);
            makeGraph(allY, allValuesYY, fieldNames[i]+"_YY", path, gFormat);
            makeGraph(allY, allValuesZZ, fieldNames[i]+"_ZZ", path, gFormat);
            makeGraph(allY, allValuesXY, fieldNames[i]+"_XY", path, gFormat);
            makeGraph(allY, allValuesXZ, fieldNames[i]+"_XZ", path, gFormat);
            makeGraph(allY, allValuesYZ, fieldNames[i]+"_YZ", path, gFormat);
        }
        else
        {
            Info << "    Unsupported field, ignoring field " << fieldNames[i] << endl;
            Info << "    Please add this field in collapse.H " << endl;
        }

    }
