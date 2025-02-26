#include "Parameters.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

Parameters::Parameters(std::string& inputFileName, std::string& inputRestartFolder, std::string& runModeOverWrite, std::string& restartStepOverwrite, std::string& outputSubFolderExplicit ) {
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
        }else if (a=="startingStepNum") {
            split >> l;
            startingStepNum = l;
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
            split >> l;
            FIRE_N_positive_min = l;
        }else if (a=="FIRE_N_negative_max") {
            split >> l;
            FIRE_N_negative_max = l;
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
        }else if (a=="FIRE_RTolerance") {
            split >> c;
            FIRE_RTolerance = c;
        }else if (a=="FIRE_numStrainSteps") {
            split >> d;
            FIRE_numStrainSteps = d;
        }else if (a=="FIRE_startingStrainStep") {
            split >> l;
            FIRE_startingStrainStep = l;
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
        }else if (a=="dumpSmoothenCurves") {
        split >> trueFalse;
        dumpSmoothenCurves = trueFalse;
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
        }else if (a=="ntnRepulsionMethod") {
        split >> f;
        ntnRepulsionMethod = f;
        }else if (a=="ntsHarmonicPenaltyStiffness") {
            split >> c;
            ntsHarmonicPenaltyStiffness = c;
        }else if (a=="ntsPowerlawRepulseEnergy") {
        split >> c;
        ntsPowerlawRepulseEnergy = c;
        }else if (a=="ntnPLEnergy") {
        split >> c;
        ntnPLEnergy = c;
        }else if (a=="ntsPowerlawLjScale") {
            split >> c;
            ntsPowerlawLjScale = c;
        }else if (a=="ntnRadius") {
        split >> c;
        ntnRadius = c;
        }else if (a=="ntnHStiffness") {
        split >> c;
        ntnHStiffness = c;
        }else if (a=="ntnPLRcutoffOverRadius") {
        split >> c;
        ntnPLRcutoffOverRadius = c;
        }else if (a=="integrator") {
        split >> d;
        integrator = d;
        }else if (a=="gntn_NGhostNodes") {
        split >> d;
        gntn_NGhostNodes = d;
        }else if (a=="targetPressure") {
            split >> c;
            targetPressure = c;
        }else if (a=="shearTo") {
            split >> c;
            shearTo = c;
        }else if (a=="writeToConsoleEvery") {
            split >> l;
            writeToConsoleEvery = l;
        }else if (a=="smoothCorners") {
            split >> trueFalse;
            smoothCorners = trueFalse;
        }else if (a=="alpha_HermitPol") {
            split >> c;
            alpha_HermitPol = c;
        }
        
        
    }
    //    ioverwtire are parameter recieved as console argumentsa
    if (runModeOverWrite=="stepShear"){
        runMode = "stepShear";
        startingMode = "restart";
        startingStepNum = std::stoi(restartStepOverwrite);
        restartFile = inputRestartFolder+"/data-per-node-"+restartStepOverwrite+".txt";
        outputSubFolder = outputSubFolderExplicit;
        
    }
    if (runModeOverWrite=="surfaceShear"){
        runMode = "surfaceShear";
        startingMode = "restart";
        startingStepNum = std::stoi(restartStepOverwrite);
        restartFile = inputRestartFolder+"/data-per-node-"+restartStepOverwrite+".txt";
        outputSubFolder = outputSubFolderExplicit;
        
    }
    if (runModeOverWrite=="contineousShear"){
        runMode = "contineousShear";
        startingMode = "restart";
        startingStepNum = std::stoi(restartStepOverwrite);
        restartFile = inputRestartFolder+"/data-per-node-"+restartStepOverwrite+".txt";
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
    printit("Ap",Ap);
    
    printit("dumpEvery",dumpEvery);
    printit("writeToConsoleEvery",writeToConsoleEvery);
    printit("splitDataEvery",splitDataEvery);
    printit("dumpSmoothenCurves",dumpSmoothenCurves);
    printit("startingTimeStep",startingStepNum);
    printit("outputFolderName",outputFolderName);
    printit("surfaceNodesFileName",surfaceNodesFileName);
    printit("trianglesFileName",trianglesFileName);
    printit("initialNodesFileName",initialNodesFileName);
    
    printit("initTopPos",initTopPos);
    printit("initBotPos",initBotPos);
    printit("initRightPos",initRightPos);
    printit("initLeftPos",initLeftPos);
    
    
    
    
    printit("solver",solver);
    if (solver=="GD"){
        printit("deformationRate",deformationRate);
        printit("dt",dt);
        printit("maxForceTol",maxForceTol);
    }else if (solver=="FIRE2" || solver=="FIRE"){
        printit("FIRE_dtmax",FIRE_dtmax);
        printit("FIRE_Nmax",FIRE_Nmax);
        printit("FIRE_finc",FIRE_finc);
        printit("FIRE_fdec",FIRE_fdec);
        printit("FIRE_alpha_start",FIRE_alpha_start);
        printit("FIRE_falpha",FIRE_falpha);
        printit("FIRE_dt_start",FIRE_dt_start);
        printit("FIRE_RTolerance",FIRE_RTolerance);
        
        if (solver=="FIRE2"){
            printit("FIRE_dtmin",FIRE_dtmin);
            printit("FIRE_N_negative_max",FIRE_N_negative_max);
            printit("FIRE_intialdelay",FIRE_intialdelay);
        }
    }
    
    printit("integrator", integrator);
    
    printit("runMode",runMode);
    if (runMode=="compress"){
        printit("targetPhi",targetPhi);
        if (startingMode=="new") printit("initialStretch",initialStretch);
        if (solver=="FIRE2" || solver=="FIRE"){
            printit("FIRE_numStrainSteps",FIRE_numStrainSteps);
            printit("FIRE_startingStrainStep",FIRE_startingStrainStep);
        }
    }else if (runMode=="stepShear" || runMode=="surfaceShear" || runMode=="continuousShear"){
        printit("targetShear",targetShear);
    }else if (runMode=="special"){
        printit("targetPressure",targetPressure);
        printit("shearTo",shearTo);
    }
    
    printit("startingMode",startingMode);
    if (startingMode=="restart"){
        printit("restartFile",restartFile);
    }
    
    printit("boundaryType",boundaryType);
    if (boundaryType=="walls"){
        printit("wallStyle",wallStyle);
        printit("HWallStiffness",HWallStiffness);
        printit("PLWallEnergyScale",PLWallEnergyScale);
        printit("PLWallLJScale",PLWallLJScale);
    }else if (boundaryType=="periodic"){
        printit("imagesMargin",imagesMargin);
        printit("dumpPeriodicImagesXY",dumpPeriodicImagesXY);
        
    }
    
    printit("contactMethod", contactMethod);
    if (contactMethod=="nts"){
        printit("ntsHarmonicPenaltyStiffness",ntsHarmonicPenaltyStiffness);
        printit("segmentCellMethod",segmentCellMethod);
        printit("identifyAndDumbFacets",identifyAndDumbFacets);
        printit("reversibleMasterSlaveRole", reversibleMasterSlaveRole);
        printit("smoothCorners", smoothCorners);
        printit("alpha_HermitPol", alpha_HermitPol);
        

    }else if (contactMethod=="ntn" || contactMethod=="gntn"){
        if (contactMethod=="gntn") printit("gntn_NGhostNodes", gntn_NGhostNodes);
        printit("ntnRepulsionMethod", ntnRepulsionMethod);
        if (ntnRepulsionMethod == "powerlaw"){
            printit("ntnPLEnergy", ntnPLEnergy);
            printit("ntnRadius", ntnRadius);
            printit("ntnPLRcutoffOverRadius", ntnPLRcutoffOverRadius);
        }else if (ntnRepulsionMethod == "harmonic"){
            printit("ntnHStiffness", ntnHStiffness);
            printit("ntnRadius", ntnRadius);
        }
    }
    printit("verletCellCutoff",verletCellCutoff);
    
    printit("calculateHessian",calculateHessian);

            std::cout << std::endl;
}


    template <typename T>
void Parameters::printit(const std::string name,const T v) const {
        std::cout << name+" = " << v << std::endl;
        std::ofstream logfile;
        logfile.open (runMode+"-parameters.log", std::ios_base::app);
        logfile << name+" = " << v << std::endl;
        logfile.close();
        
        std::ofstream logfile2;
        logfile2.open (outputFolderName+"/used-parameters.log", std::ios_base::app);
        logfile2 << name+" = " << v << std::endl;
        logfile2.close();
    }

