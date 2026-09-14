    // Note: if you want to add a new addModel, do it here

    scalar thrSurf(readScalar(HFDIBDEMDict_.lookup("surfaceThreshold")));

    forAll (bodyNames_,nameI)
    {
        word bodyName(bodyNames_[nameI]);
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "Creating addition model for body: " << bodyName << endl;

        // initialize the addition models
        dictionary addModelDict(HFDIBDEMDict_.subDict(bodyName).subDict("bodyAddition"));
        word howToAdd(word(addModelDict.lookup("addModel")));
        word bodyGeom;
        std::unique_ptr<geomModel> bodyGeomModel;
        // check if the immersedDict_ contains bodyGeom
        if (HFDIBDEMDict_.subDict(bodyName).found("bodyGeom"))
        {
            const word input(HFDIBDEMDict_.subDict(bodyName).lookup("bodyGeom"));
            bodyGeom = input;
            InfoH << iB_Info << "Found bodyGeom for " << bodyName
                << ", the body is: " << bodyGeom << endl;
        }
        else
        {
            notImplemented("stop");
            bodyGeom = "convex";
            InfoH << iB_Info << "Did not find bodyGeom for "
                << bodyName << ", using bodyGeom: " << bodyGeom << endl;
        }

        {
            word stlPath("constant/triSurface/"+bodyName+".stl");
            InfoH << iB_Info << "bodyGeom: " << bodyGeom
                << " not supported, using bodyGeom nonConvex" << endl;
            bodyGeom = "nonConvex";
            bodyGeomModel = std::unique_ptr<nonConvexBody>(
                 new nonConvexBody(mesh_,stlPath,thrSurf));
        }

        InfoH << addModel_Info << "-- addModelMessage-- "
            << "trying to set addModel: " << howToAdd << endl;

        if (howToAdd == "once")
        {
            addModels_.set
            (
                nameI,
                new addModelOnce
                (
                    addModelDict,
                    mesh_,
                    startTime0,
                    std::move(bodyGeomModel),
                    cellPoints_
                )
            );
        }
        else
        {
            notImplemented("STOP");
        }
    }
