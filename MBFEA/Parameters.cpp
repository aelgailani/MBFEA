#include "Parameters.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

Parameters::Parameters(std::string& inputFileName, std::string& inputRestartFolder, std::string& runModeOverWrite, std::string& restartStepOverwrite) {
    std::ifstream inFile;
    inFile.open(inputFileName); //Make sure later that the file is open in "read" mode only
    //Check for errors
    if (inFile.fail())
    {
        std::cerr << "Error Openning Input Parameters File" << std::endl;
        exit(1);
    }
    
    std::string line;
    std::string a,b,f;
    double c;
    int d;
    long l;
    bool trueFalse;
    while (std::getline(inFile, line)) {
        std::istringstream split(line);
        split >> a;
        if (a=="kTOverOmega") {
            split >> c;
            kTOverOmega = c;
        }else if (a=="NkT") {
            split >> c;
            NkT = c;
        }else if (a=="chi") {
            split >> c;
            chi = c;
        }else if (a=="verletCellCutoff") {
            split >> c;
            verletCellCutoff = c;
        }else if (a=="initialStretch") {
            split >> c;
            initialStretch = c;
        }else if (a=="dumpEvery") {
            split >> l;
            dumpEvery = l;
        }else if (a=="splitDataEvery") {
            split >> l;
            splitDataEvery = l;
        }else if (a=="startingTimeStep") {
            split >> l;
            startingTimeStep = l;
        }else if (a=="outputFolderName") {
            split >> f;
            outputFolderName = f;
        }else if (a=="surfaceNodesFileName") {
            split >> b;
            surfaceNodesFileName = b;
        }else if (a=="trianglesFileName") {
            split >> b;
            trianglesFileName = b;
        }else if (a=="initialNodesFileName") {
            split >> b;
            initialNodesFileName = b;
        }else if (a=="boundaryType") {
            split >> b;
            boundaryType = b;
        }else if (a=="imagesMargin") {
            split >> c;
            imagesMargin = c;
        }else if (a=="initTopPos") {
            split >> c;
            initTopPos = c;
        }else if (a=="initBotPos") {
            split >> c;
            initBotPos = c;
        }else if (a=="initRightPos") {
            split >> c;
            initRightPos = c;
        }else if (a=="initLeftPos") {
            split >> c;
            initLeftPos = c;
        }else if (a=="deformationRate") {
            split >> c;
            deformationRate = c;
        }else if (a=="targetPhi") {
            split >> c;
            targetPhi = c;
        }else if (a=="dt") {
            split >> c;
            dt = c;
        }else if (a=="wallStyle") {
            split >> b;
            wallStyle = b;
        }else if (a=="HWallStiffness") {
            split >> c;
            HWallStiffness = c;
        }else if (a=="PLWallEnergyScale") {
            split >> c;
            PLWallEnergyScale = c;
        }else if (a=="PLWallLJScale") {
            split >> c;
            PLWallLJScale = c;
        }else if (a=="Ap") {
            split >> c;
            Ap = c;
        }else if (a=="runMode") {
            split >> b;
            runMode = b;
        }else if (a=="startingMode") {
            split >> b;
            startingMode = b;
        }else if (a=="restartFile") {
            split >> b;
            restartFile = b;
        }else if (a=="targetShear") {
            split >> c;
            targetShear = c;
        }else if (a=="maxForceTol") {
            split >> c;
            maxForceTol = c;
        }
        else if (a=="solver") {
            split >> b;
            solver = b;
        }
        else if (a=="FIRE_dtmax") {
            split >> c;
            FIRE_dtmax = c;
        }
        else if (a=="FIRE_N_positive_min") {
            split >> c;
            FIRE_N_positive_min = c;
        }else if (a=="FIRE_N_negative_max") {
            split >> c;
            FIRE_N_negative_max = c;
        }else if (a=="FIRE_Nmax") {
            split >> l;
            FIRE_Nmax = l;
        }
        else if (a=="FIRE_finc") {
            split >> c;
            FIRE_finc = c;
        }
        else if (a=="FIRE_fdec") {
            split >> c;
            FIRE_fdec = c;
        }else if (a=="FIRE_alpha_start") {
            split >> c;
            FIRE_alpha_start = c;
        }else if (a=="FIRE_falpha") {
            split >> c;
            FIRE_falpha = c;
        }else if (a=="FIRE_dt_start") {
            split >> c;
            FIRE_dt_start = c;
        }else if (a=="FIRE_dtmin") {
            split >> c;
            FIRE_dtmin = c;
        }else if (a=="FIRE_intialdelay") {
        split >> trueFalse;
        FIRE_intialdelay = trueFalse;
        }else if (a=="RTolerance") {
            split >> c;
            RTolerance = c;
        }else if (a=="numStrainSteps") {
            split >> d;
            numStrainSteps = d;
        }else if (a=="startingStrainStep") {
            split >> l;
            startingStrainStep = l;
        }else if (a=="targetNodes") {
            while(split >> d){
                targetNodes.push_back(d);
//                std::cout << d <<"\t";
            }
        }else if (a=="gammaX") {
            split >> c;
            gammaX = c;
        }else if (a=="gammaY") {
            split >> c;
            gammaY = c;
        }else if (a=="segmentCellMethod") {
            split >> d;;
            segmentCellMethod = d;
        }else if (a=="dumpPeriodicImagesXY") {
        split >> trueFalse;
        dumpPeriodicImagesXY = trueFalse;
        }else if (a=="callPythonPlot") {
        split >> trueFalse;
        callPythonPlot = trueFalse;
        }else if (a=="calculateHessian") {
        split >> trueFalse;
        calculateHessian = trueFalse;
        }else if (a=="identifyAndDumbFacets") {
        split >> trueFalse;
        identifyAndDumbFacets = trueFalse;
        }else if (a=="reversibleMasterSlaveRole") {
        split >> trueFalse;
        reversibleMasterSlaveRole = trueFalse;
        }else if (a=="contactMethod") {
        split >> f;
        contactMethod = f;
        }else if (a=="ntsPenaltyMethod") {
        split >> f;
        ntsPenaltyMethod = f;
        }else if (a=="ntsHarmonicPenaltyStiffness") {
            split >> c;
            ntsHarmonicPenaltyStiffness = c;
        }else if (a=="ntsPowerlawRepulseEnergy") {
        split >> c;
        ntsPowerlawRepulseEnergy = c;
        }else if (a=="ntnRepulseEnergy") {
        split >> c;
        ntnRepulseEnergy = c;
        }else if (a=="ntsPowerlawLjScale") {
            split >> c;
            ntsPowerlawLjScale = c;
        }else if (a=="ntnLjScale") {
        split >> c;
        ntnLjScale = c;
        }
        
        
    }
    //    identifyAndDumbFacets
    if (runModeOverWrite=="stepShear"){
        runMode = "stepShear";
        startingMode = "restart";
        startingTimeStep = std::stoi(restartStepOverwrite);
        restartFile = inputRestartFolder+"/data_per_node-"+restartStepOverwrite+".txt";
        
    }
    if (runModeOverWrite=="contineousShear"){
        runMode = "contineousShear";
        startingMode = "restart";
        startingTimeStep = std::stoi(restartStepOverwrite);
        restartFile = inputRestartFolder+"/data_per_node-"+restartStepOverwrite+".txt";
    }
    if (outputFolderName=="auto"){
        outputFolderName = runMode;
    }
    if (outputFolderName==inputFileName){
        std::cout << "use different output directory name to avoid overwriting your data ! " << std::endl;;
            exit(1);
           }
    if (boundaryType=="walls"){
        imagesMargin = 0.0 ;
    }
}



void Parameters::print_to_console(void) const {
    printit("kTOverOmega", kTOverOmega);
    printit("NkT",NkT);
    printit("chi",chi);
    printit("verletCellCutoff",verletCellCutoff);
    printit("initialStretch",initialStretch);
    printit("dumpEvery",dumpEvery);
    printit("splitDataEvery",splitDataEvery);
    printit("startingTimeStep",startingTimeStep);
    printit("outputFolderName",outputFolderName);
    printit("surfaceNodesFileName",surfaceNodesFileName);
    printit("trianglesFileName",trianglesFileName);
    printit("initialNodesFileName",initialNodesFileName);
    printit("boundaryType",boundaryType);
    printit("imagesMargin",imagesMargin);
    printit("initTopPos",initTopPos);
    printit("initBotPos",initBotPos);
    printit("initRightPos",initRightPos);
    printit("initLeftPos",initLeftPos);
    printit("deformationRate",deformationRate);
    printit("targetPhi",targetPhi);
    printit("dt",dt);
    printit("wallStyle",wallStyle);
    printit("HWallStiffness",HWallStiffness);
    printit("PLWallEnergyScale",PLWallEnergyScale);
    printit("PLWallLJScale",PLWallLJScale);
    printit("penaltyStiffness",ntsHarmonicPenaltyStiffness);
    printit("Ap",Ap);
    printit("runMode",runMode);
    printit("startingMode",startingMode);
    printit("restartFile",restartFile);
    printit("targetShear",targetShear);
    printit("maxForceTol",maxForceTol);
    printit("solver",solver);
    printit("numStrainSteps",numStrainSteps);
    printit("startingStrainStep",startingStrainStep);
    printit("segmentCellMethod",segmentCellMethod);
    printit("dumpPeriodicImagesXY",dumpPeriodicImagesXY);
    printit("callPythonPlot",callPythonPlot);
    printit("calculateHessian",callPythonPlot);
    printit("identifyAndDumbFacets",identifyAndDumbFacets);
    printit("reversibleMasterSlaveRole", reversibleMasterSlaveRole);
    printit("contactMethod", contactMethod);
    printit("repulseEnergy", ntnRepulseEnergy);
    printit("ljScale", ntnLjScale);


            std::cout << std::endl;
}


    template <typename T>
void Parameters::printit(const std::string name,const T v) const {
        std::cout << name+" = " << v << std::endl;
}
