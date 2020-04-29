#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <cmath>
#include <sys/stat.h>
#include <dirent.h>
#include <iomanip>
#include <valarray>
#include "Parameters.hpp"
#include "Configuration.hpp"
#include "BaseSysData.hpp"

Configuration::Configuration(const BaseSysData& baseData, const Parameters& pars){
    
    if (pars.startingMode=="new") {
        
        curPosX = baseData.refPosX * pars.initialStretch;
        curPosY = baseData.refPosY * pars.initialStretch;
        topPos = pars.initTopPos;
        botPos = pars.initBotPos;
        rightPos = pars.initRightPos;
        leftPos = pars.initLeftPos;
        lyCur = baseData.lyRef;
        lxCur =  baseData.lxRef;
        lyNew = lyCur;
        lxNew = lxCur;


    }else if (pars.startingMode=="restart")  {
        std::ifstream inFile;
        inFile.open(pars.restartFile); //Make sure later that the file is open in "read" mode only
        //Check for errors
        if (inFile.fail())
        {
            std::cerr << "Error Openning Restart Nodes File" << std::endl;
            exit(1);
        }
        std::string line,a;
        double b,c,d,x,y;
        long l;
        std::getline(inFile, line);
        std::istringstream split1(line);
        split1 >> a;
        
        if (a == "Basic_data:") {
        
            while (a != "Nodes_data:"){
                std::getline(inFile, line);
                std::istringstream split(line);
                split >> a;
        
                if (a=="numNodes") {
                    split >> l;
                    assert(l== baseData.numOriginalNodes);
                }
            
                if (a=="kTOverOmega") {
                    split >> b;
                    assert(b==pars.kTOverOmega);
                }
                
                if (a=="timeStep") {
                    split >> l;
                }

                if (a=="original_box_LRBT") {
                    split >> b >> c >> x >> y;
                    leftPos = b;
                    rightPos = c;
                    botPos = x;
                    topPos = y;
                }

            }
            
            lxCur = rightPos - leftPos;
            lyCur = topPos - botPos;
            lyNew = lyCur;
            lxNew = lxCur;
            std::getline(inFile, line);
            int id = 0;
            curPosX.resize(baseData.numOriginalNodes);
            curPosY.resize(baseData.numOriginalNodes);
            while(std::getline(inFile, line))
            {
                std::istringstream split(line);
                split >> b >> x >> y >> c >> d;
                curPosX(id) = x;
                curPosY(id) = y;
                id++;
            }
            assert(baseData.numOriginalNodes==curPosX.size());
            assert(baseData.numOriginalNodes==curPosY.size());
            inFile.close();
        
        //if the restart file is in the old format, you can remove this later option later.
        }else{
                
                if (a=="numNodes") {
                    split1 >> b;
                    assert(b== baseData.numOriginalNodes);
                }
                
                std::getline(inFile, line);
                std::istringstream split2(line);
                split2 >> a;
                std::cout << a << std::endl;
                if (a=="kTOverOmega") {
                    split2 >> b;
                    assert(b==pars.kTOverOmega);
                    
                }
                
                std::getline(inFile, line); //this reads phi but skip it
                std::getline(inFile, line);
                std::istringstream split3(line);
                split3 >> a;
                if (a=="timeStep") {
                    split3 >> b;
                }
                std::getline(inFile, line);
                std::istringstream split4(line);
                split4 >> a;
                std::cout << a << std::endl;
                if (a=="wallsLRBT") {
                    split4 >> b >> c >> x >> y;
                    leftPos = b;
                    rightPos = c;
                    botPos = x;
                    topPos = y;
                }
                
                lxCur = rightPos - leftPos;
                lyCur = topPos - botPos;
                lyNew = lyCur;
                lxNew = lxCur;
            
                std::getline(inFile, line);
                std::getline(inFile, line);
                int id = 0;
                
                curPosX.resize(baseData.numOriginalNodes);
                curPosY.resize(baseData.numOriginalNodes);
                while(std::getline(inFile, line))
                {
                    std::istringstream split(line);
                    split >> b >> x >> y >> c >> d;
                    curPosX(id) = x;
                    curPosY(id) = y;
                    id++;
                }
                assert(baseData.numOriginalNodes==curPosX.size());
                assert(baseData.numOriginalNodes==curPosY.size());
                inFile.close();
            }
        }

    ex = log(lxNew/baseData.lxRef);
    ey = log(lyNew/baseData.lyRef);
    A = lxNew * lyNew;
    e0 = - 0.5*(ex+ey); // my conviension is positive for compression
    e1 = (ex-ey);
    phi = pars.Ap / A;
    
    curPosXAtLastStep = curPosX;
    curPosYAtLastStep = curPosY;
    displacementSinceLastStep = curPosX - curPosX ;  // basically a zero vector, initially
    defGradXX.resize(baseData.numElements,1);
    defGradXY.resize(baseData.numElements,1);
    defGradYX.resize(baseData.numElements,1);
    defGradYY.resize(baseData.numElements,1);
    invDefGradTransXX.resize(baseData.numElements,1);
    invDefGradTransXY.resize(baseData.numElements,1);
    invDefGradTransYX.resize(baseData.numElements,1);
    invDefGradTransYY.resize(baseData.numElements,1);
    PK1stressXX.resize(baseData.numElements,1);
    PK1stressXY.resize(baseData.numElements,1);
    PK1stressYX.resize(baseData.numElements,1);
    PK1stressYY.resize(baseData.numElements,1);
    CstressXX.resize(baseData.numElements,1);
    CstressXY.resize(baseData.numElements,1);
    CstressYX.resize(baseData.numElements,1);
    CstressYY.resize(baseData.numElements,1);
    swellingPressurePerEle.resize(baseData.numElements,1);
    fSquared.resize(baseData.numElements,1);
    elasticEnergyPerEle.resize(baseData.numElements,1);
    mixingEnergyPerEle.resize(baseData.numElements,1);
    refArea.resize(baseData.numElements,1);
    areaRatio.resize(baseData.numElements,1);
    
    gradX.resize(baseData.numElements,baseData.numOriginalNodes);
    gradY.resize(baseData.numElements,baseData.numOriginalNodes);
    
    for (int triID=0; triID<baseData.numElements; triID++) {
        double area = 0.0;
        int a = baseData.triangles.at(triID)[0];
        int b = baseData.triangles.at(triID)[1];
        int c = baseData.triangles.at(triID)[2];
        area += baseData.refPosX[a]*baseData.refPosY[b];
        area += baseData.refPosX[b]*baseData.refPosY[c];
        area += baseData.refPosX[c]*baseData.refPosY[a];
        area -= baseData.refPosX[b]*baseData.refPosY[a];
        area -= baseData.refPosX[c]*baseData.refPosY[b];
        area -= baseData.refPosX[a]*baseData.refPosY[c];
        refArea[triID]=fabs(area)*0.5;
        gradY.insert(triID,a)=-(baseData.refPosX[b]-baseData.refPosX[c])/area;
        gradY.insert(triID,b)=-(baseData.refPosX[c]-baseData.refPosX[a])/area;
        gradY.insert(triID,c)=-(baseData.refPosX[a]-baseData.refPosX[b])/area;
        gradX.insert(triID,a)=(baseData.refPosY[b]-baseData.refPosY[c])/area;
        gradX.insert(triID,b)=(baseData.refPosY[c]-baseData.refPosY[a])/area;
        gradX.insert(triID,c)=(baseData.refPosY[a]-baseData.refPosY[b])/area;
        
    }
    
    // These are for debugging mainly
//    consistencyFactorX.resize(baseData.numOriginalNodes);
//    consistencyFactorY.resize(baseData.numOriginalNodes);
//    consistencyErrorFactorY.resize(baseData.numOriginalNodes);
//    consistencyErrorFactorX.resize(baseData.numOriginalNodes);
    //Hessian debugging
//    FXref.resize(baseData.numOriginalNodes);
//    FYref.resize(baseData.numOriginalNodes);
//    DeltaForceXRatio.resize(baseData.numOriginalNodes);
//    DeltaForceYRatio.resize(baseData.numOriginalNodes);
//    refX.resize(baseData.numOriginalNodes);
//    refY.resize(baseData.numOriginalNodes);
//    HessianFX.resize(baseData.numOriginalNodes);
//    HessianFY.resize(baseData.numOriginalNodes);
    
    // These variables are for the linearization of the second energy derivitive:
//    Wm_primePrime.resize(baseData.numElements,1);
     if (pars.calculateHessian){
    KMxxjx.resize(baseData.numElements,baseData.numOriginalNodes);
    KMxxjy.resize(baseData.numElements,baseData.numOriginalNodes);
    KMxyjx.resize(baseData.numElements,baseData.numOriginalNodes);
    KMxyjy.resize(baseData.numElements,baseData.numOriginalNodes);
    KMyxjx.resize(baseData.numElements,baseData.numOriginalNodes);
    KMyxjy.resize(baseData.numElements,baseData.numOriginalNodes);
    KMyyjx.resize(baseData.numElements,baseData.numOriginalNodes);
    KMyyjy.resize(baseData.numElements,baseData.numOriginalNodes);
    Hixjx.resize(baseData.numOriginalNodes,baseData.numOriginalNodes);
    Hixjy.resize(baseData.numOriginalNodes,baseData.numOriginalNodes);
    Hiyjx.resize(baseData.numOriginalNodes,baseData.numOriginalNodes);
    Hiyjy.resize(baseData.numOriginalNodes,baseData.numOriginalNodes);
    InvHixjx.resize(baseData.numOriginalNodes,baseData.numOriginalNodes);
    InvHixjy.resize(baseData.numOriginalNodes,baseData.numOriginalNodes);
    InvHiyjx.resize(baseData.numOriginalNodes,baseData.numOriginalNodes);
    InvHiyjy.resize(baseData.numOriginalNodes,baseData.numOriginalNodes);
    affineForceX.resize(baseData.numOriginalNodes);
    affineForceY.resize(baseData.numOriginalNodes);
    nonAffineVX.resize(baseData.numOriginalNodes);
    nonAffineVY.resize(baseData.numOriginalNodes);
    Hessian.resize(baseData.numOriginalNodes*2,baseData.numOriginalNodes*2);
    augmentedNonAffineV.resize(baseData.numOriginalNodes*2);
    augmentedAffineF.resize(baseData.numOriginalNodes*2);

     }
    
    // If slover is FIRE, initiate zero velocity vectors
    if (pars.solver == "FIRE" || pars.solver == "FIRE2"){
        velocityX.resize(baseData.numOriginalNodes);
        velocityY.resize(baseData.numOriginalNodes);
        prevVelocityX.resize(baseData.numOriginalNodes);
        prevVelocityY.resize(baseData.numOriginalNodes);
        
    }
}



void Configuration::update_post_processing_data(const BaseSysData& baseData, const Parameters& pars){
    
    maxR = sqrt((forceX.array()*forceX.array()+forceY.array()*forceY.array()).maxCoeff());
    avgR = sqrt((forceX.array()*forceX.array()+forceY.array()*forceY.array()).mean());
//    LX = (curPosX.maxCoeff()-curPosX.minCoeff());
//    LY = (curPosY.maxCoeff()-curPosY.minCoeff());
    LX = (rightPos-leftPos);
    LY = (topPos-botPos);
    Fh = 0.5*(-wallForceRight.sum()+wallForceLeft.sum());
    Fv = 0.5*(-wallForceTop.sum()+wallForceBottom.sum());
    P1 = 0.5*(Fv/lxNew+Fh/lyNew);
//    P2 = ((CstressYY+CstressYX).dot((refArea.array()*areaRatio.array()).matrix()) + (CstressXX+CstressXY).dot((refArea.array()*areaRatio.array()).matrix()))/(2*LX*LY);
    P2 = ((CstressYY).dot((refArea.array()*areaRatio.array()).matrix()) + (CstressXX).dot((refArea.array()*areaRatio.array()).matrix())+(KWoodYY+KWoodXX))/(2*LX*LY);
    
    S1 = 0.5*(Fv/lxNew-Fh/lyNew);
    S2 = ((CstressYY).dot((refArea.array()*areaRatio.array()).matrix())- (CstressXX).dot((refArea.array()*areaRatio.array()).matrix())+(KWoodYY-KWoodXX))/(2*LX*LY);
//    S2 = ((CstressYY+CstressYX).dot((refArea.array()*areaRatio.array()).matrix())- (CstressXX+CstressXY).dot((refArea.array()*areaRatio.array()).matrix()))/(2*LX*LY);
    ex = log(lxNew/baseData.lxRef);
    ey = log(lyNew/baseData.lyRef);
    A = lxNew * lyNew;
    A_material = refArea.dot(areaRatio);
    e0 = - 0.5*(ex+ey); // my conviension is positive for compression
    e1 = (ex-ey);
    phi = pars.Ap / A;
    deltaTotEnergy = totalEnergy - prevTotEnergy;
    prevTotEnergy = totalEnergy;
    L2NormResidual = - (forceX.dot(forceX) + forceY.dot(forceY));
    
    shearVirial = ((CstressYY+CstressYX).dot((refArea.array()*areaRatio.array()).matrix())- (CstressXX+CstressXY).dot((refArea.array()*areaRatio.array()).matrix()));
    pressureVirial = ((CstressYY+CstressYX).dot((refArea.array()*areaRatio.array()).matrix())+ (CstressXX+CstressXY).dot((refArea.array()*areaRatio.array()).matrix()));
    
    DPOverDe0 = (prev_P2-P2)/(prev_e0-e0);
    DSOverDe1 = (prev_S2-S2)/(prev_e1-e1);
    prev_e0 = e0;
    prev_e1 =e1;
    prev_P2 = P2;
    prev_S2 = S2;
}

void Configuration::dump_global_data(const Parameters& pars, const long& timeStep, std::string mode, std::string purpose){
    std::string fname ;
    std::string first = std::to_string(timeStep/pars.splitDataEvery*pars.splitDataEvery);
    std::string last =std::to_string(timeStep/pars.splitDataEvery*pars.splitDataEvery+pars.splitDataEvery-1);
    
    
    if (purpose=="running"){
        fname = "/data-steps-"+first+"-"+last+".txt";
    }else if (purpose=="final"){
        fname = "/final-data.txt";
    }
   
    if (mode=="write" || first != lastStepFirst){
        std::ofstream dataFile;
        if (pars.runMode == "stepShear" || pars.runMode == "continuousShear"){
        dataFile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingTimeStep))).c_str()+fname);
        }else{
            dataFile.open (pars.outputFolderName+fname);
        }
        dataFile
        <<  "step"  << std::setw(30)
        <<  "totalEnergy"  << std::setw(20)
        <<  "contactEnergy"  << std::setw(20)
        <<  "shearVirial"  << std::setw(20)
        <<  "pressureVirial"  << std::setw(20)
        <<  "maxResidualF"  << std::setw(20)
        <<  "pressure1"  << std::setw(20)
        <<  "pressure2"  << std::setw(20)
        <<  "KWpressure2"  << std::setw(20)
        <<  "shearStress1"  << std::setw(20)
        <<  "shearStress2"  << std::setw(20)
        <<  "KWshearStress2"  << std::setw(20)
        <<  "boxArea"  << std::setw(20)
        <<  "materialArea"  << std::setw(20)
        <<  "phi"  << std::setw(20)
        <<  "e0"  << std::setw(20)
        <<  "e1" << std::setw(20)
        <<  "DPOverDe0" << std::setw(20)
        <<  "DSOverDe1qq" << std::setw(20)
        <<  "dt" << std::setw(20)
        <<  "defRate" << std::setw(20)
        <<  "penalty" << std::endl;
        dataFile.close();
        
        
    }else if (mode=="append"){
        std::ofstream dataFile;
        if (pars.runMode == "stepShear" || pars.runMode == "continuousShear"){
        dataFile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingTimeStep))).c_str()+fname,  std::ios_base::app);
        }else{
            dataFile.open (pars.outputFolderName+fname,  std::ios_base::app);
        }
        dataFile << std::setprecision(9)
        <<  timeStep  << std::setw(30)
        <<  totalEnergy  << std::setw(20)
        <<  contactsEnergy  << std::setw(20)
        <<  shearVirial  << std::setw(20)
        <<  pressureVirial  << std::setw(20)
        <<  maxR  << std::setw(20)
        <<  P1  << std::setw(20)
        <<  P2  << std::setw(20)
        <<  (KWoodYY+KWoodXX)/(2*LX*LY)  << std::setw(20)
        <<  S1  << std::setw(20)
        <<  S2  << std::setw(20)
        <<  (KWoodYY-KWoodXX)/(2*LX*LY)   << std::setw(20)
        <<  A  << std::setw(20)
        <<  A_material  << std::setw(20)
        <<  phi  << std::setw(20)
        <<  e0  << std::setw(20)
        <<  e1 << std::setw(20)
        <<  DPOverDe0 << std::setw(20)
        <<  DSOverDe1 << std::setw(20)
        <<  pars.dt  << std::setw(20)
        <<  pars.deformationRate  << std::setw(20)
        <<  pars.ntsHarmonicPenaltyStiffness << std::endl;
        dataFile.close();
        
        
    }
    
   lastStepFirst = first;
}


void Configuration::dump_per_node(const BaseSysData& baseData, const Parameters& pars, long& timeStep){
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    if (pars.runMode == "stepShear" || pars.runMode == "continuousShear"){
    myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingTimeStep))+"/data-per-node-"+step+".txt").c_str());
    }else{
        myfile.open (pars.outputFolderName+"/data-per-node-"+step+".txt");
    }
    
    myfile << "Basic_data:" << std::endl;
    myfile << "timeStep" << "\t" << timeStep << std::endl;
    myfile << "number_of_nodes" << "\t" << baseData.numOriginalNodes << std::endl;
    myfile << "kTOverOmega" << "\t" << pars.kTOverOmega << std::endl;
    myfile << "phi" << "\t" << phi << std::endl;
    myfile << "e0_strain" << "\t" << e0 << std::endl;
    myfile << "e1_strain" << "\t" << e1 << std::endl;
    myfile << "dt" << "\t" << pars.dt << std::endl;
    myfile << "max_residual" << "\t" << maxR << std::endl;
    myfile << "deformation_rate" << "\t" << pars.deformationRate   << std::endl;
    myfile << "penalty_stiffness" << "\t" << pars.ntsHarmonicPenaltyStiffness << std::endl;
    myfile << "verlet_cell_cutoff" << "\t" << pars.verletCellCutoff << std::endl;
    myfile << "original_box_LRBT" << "\t" << leftPos << "\t" << rightPos << "\t" << botPos << "\t" << topPos<< "\n"  << std::endl;
    myfile << "Nodes_data:" << std::endl;
    myfile
    << "id"  << std::setw(20)
    <<  "x"  << std::setw(20)
    <<  "y"  << std::setw(20)
    <<  "fx"  << std::setw(20)
    <<  "fy"  << std::setw(20)
    <<  "contactFx"  << std::setw(20)
    <<  "contactFy"  << std::endl;
    
    for (int i=0; i<curPosX.size();i++ ) {
        myfile
        <<  std::setprecision(9)
        << i << std::setw(20)
        << curPosX[i] << std::setw(20)
        << curPosY[i] << std::setw(20)
        << forceX[i] << std::setw(20)
        << forceY[i] << std::setw(20)
        << surfaceForceX[i] << std::setw(20)
        << surfaceForceY[i] << std::endl;
    }
    
    myfile.close();

}

void Configuration::dump_per_node_periodic_images_on(const BaseSysData& baseData, const Parameters& pars, long& timeStep){
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    if (pars.runMode == "stepShear" || pars.runMode == "continuousShear"){
    myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingTimeStep))+"/data-per-node-periodic-"+step+".txt").c_str());
    }else{
        myfile.open (pars.outputFolderName+"/data-per-node-periodic-"+step+".txt");
    }
    myfile << "Basic_data:" << std::endl;
    myfile << "timeStep" << "\t" << timeStep << std::endl;
    myfile << "number_of_augmented_nodes" << "\t" << baseData.numNodes << std::endl;
    myfile << "images_margin" << "\t" << pars.imagesMargin << std::endl;
    myfile << "effective_simulation_box_left_boundary" << "\t" << leftPos - pars.imagesMargin*lxNew << std::endl;
    myfile << "effective_simulation_box_bot_boundary" << "\t" << botPos - pars.imagesMargin*lyNew << std::endl;
    myfile << "numCellsX" << "\t" << numXCells << std::endl;
    myfile << "numCellsY" << "\t" << numYCells << std::endl;
    myfile << "cellSizeX" << "\t" << verletCellSizeX << std::endl;
    myfile << "cellSizeY" << "\t" << verletCellSizeY <<"\n" << std::endl;
   
    myfile << "Nodes_data:" << std::endl;
    myfile
    << "id"  << std::setw(20)
    <<  "x"  << std::setw(20)
    <<  "y"  << std::setw(20)
    <<  "fx"  << std::setw(20)
    <<  "fy"  << std::setw(20)
    <<  "contactFx"  << std::setw(20)
    <<  "contactFy"  << std::endl;
    
    for (int i=0; i<augmentedCurPosX.size();i++ ) {
        myfile
        <<  std::setprecision(9)
        << i << std::setw(20)
        << augmentedCurPosX[i] << std::setw(20)
        << augmentedCurPosY[i] << std::setw(20)
        << forceX[i % baseData.numOriginalNodes] << std::setw(20)
        << forceY[i % baseData.numOriginalNodes] << std::setw(20)
        << surfaceForceX[i % baseData.numOriginalNodes] << std::setw(20)
        << surfaceForceY[i % baseData.numOriginalNodes] << std::endl;
    }
    myfile.close();
    
}


void Configuration::dump_per_ele(const BaseSysData& baseData, const Parameters& pars, long& timeStep){
    
    
    // calculated strains
    DVxDx = gradX * forceX;
    DVxDy = gradY * forceX;
    DVyDx = gradY * forceY;
    DVyDy = gradY * forceY;
    
    
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    if (pars.runMode == "stepShear" || pars.runMode == "continuousShear"){
    myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingTimeStep))+"/data-per-ele-"+step+".txt").c_str());
    }else{
        myfile.open (pars.outputFolderName+"/data-per-ele-"+step+".txt");
    }
    myfile << "Basic_data:" << std::endl;
    myfile << "timeStep" << "\t" << timeStep << std::endl;
    myfile << "number_of_elements" << "\t" << baseData.numElements << std::endl;
    myfile << "internal_energy" << "\t" <<  std::setprecision(9)<< internalEnergy << std::endl;
    myfile << "contacts_energy" << "\t" <<  std::setprecision(9)<< contactsEnergy << std::endl;
    myfile << "total_energy" << "\t" <<  std::setprecision(9)<< totalEnergy <<"\n"<< std::endl;
    myfile << "Elements_data:" << std::endl;
    myfile
    <<"id" << std::setw(20)
    << "refArea" << std::setw(20)
    << "areaRatio" << std::setw(20)
    << "FXX" << std::setw(20)
    << "FXY" << std::setw(20)
    << "FYX" << std::setw(20)
    << "FYY" << std::setw(20)
    << "PK1StressXX" << std::setw(20)
    << "PK1StressXY" << std::setw(20)
    << "PK1StressYX" << std::setw(20)
    << "PK1StressYY" << std::setw(20)
    << "CStressXX" << std::setw(20)
    << "CStressXY" << std::setw(20)
    << "CStressYX" << std::setw(20)
    << "CStressYY" << std::setw(20)
    << "DVxDx" << std::setw(20)
    << "DVxDy" << std::setw(20)
    << "DVyDx" << std::setw(20)
    << "DVyDy" << std::setw(20)
    << "swellingPressure" << std::setw(20)
    << "elasticEnergy" << std::setw(20)
    << "mixingEnergy" << std::setw(20)
    << "internalEnergy" << std::endl;

    for (int i=0; i<refArea.size();i++ ) {
        myfile
        << i << std::setw(20)
        << refArea[i] << std::setw(20)
        << areaRatio[i] << std::setw(20)
        << defGradXX[i] << std::setw(20)
        << defGradXY[i] << std::setw(20)
        << defGradYX[i] << std::setw(20)
        << defGradYY[i] << std::setw(20)
        << PK1stressXX[i] << std::setw(20)
        << PK1stressXY[i] << std::setw(20)
        << PK1stressYX[i] << std::setw(20)
        << PK1stressYY[i] << std::setw(20)
        << CstressXX[i] << std::setw(20)
        << CstressXY[i] << std::setw(20)
        << CstressYX[i] << std::setw(20)
        << CstressYY[i] << std::setw(20)
        << DVxDx[i] << std::setw(20)
        << DVxDy[i] << std::setw(20)
        << DVyDx[i] << std::setw(20)
        << DVyDy[i] << std::setw(20)
        << swellingPressurePerEle[i] << std::setw(20)
        << elasticEnergyPerEle[i] << std::setw(20)
        << mixingEnergyPerEle[i] << std::setw(20)
        << internalEnergyPerEle[i] << std::endl;
    }
    
    myfile.close();

}



void Configuration::dump_facets(const BaseSysData& baseData, const Parameters& pars, long& timeStep){
    
    
    //    ifstream stream(file);
    //    for(auto& kv : stored) {
    //      stream << kv.second << '\n';
    //      // Add '\n' character  ^^^^
    //    }
    //    stream.close();
    
    
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    if (pars.runMode == "stepShear" || pars.runMode == "continuousShear"){
    myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingTimeStep))+"/facets-"+step+".txt").c_str());
    }else{
        myfile.open (pars.outputFolderName+"/facets-"+step+".txt");
    }
    myfile << "Basic_data:" << std::endl;
    myfile << "timeStep" << "\t" << timeStep << std::endl;
    myfile << "\n" << std::endl;
    myfile << "Facets_data: (in the smNodes, the first entery is a slave node followed by its masters or master nodes or node, then the next slave node and so forth. Notice that slave nodes never duplicate in a row but master nodes can be shared and hence duplicated. Each slave node has two master nodes at max at one master mesh and a minimum of one.)" << std::endl;
    myfile
    <<"sMesh" << std::setw(7)
    << "mMesh" << std::setw(7)
    << "smNodes" <<  std::endl;

    for(auto& key : facets) {
        myfile << std::to_string(key.first.first) << std::setw(7) << std::to_string(key.first.second) << std::setw(7);
        
        for (auto& node: key.second) {
            myfile << std::to_string(node) << std::setw(7);
        }
        myfile << std::setw(-7)<< std::endl;
    }
    myfile << "EOF";
    myfile.close();

}



void Configuration::affine_compression(const BaseSysData& baseData, const Parameters& pars, double strain){
    
    std::cout << "**** compressing ****   rate " << pars.deformationRate << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
    
    yMid=0.5*(topPos+ botPos);
    xMid=0.5*(rightPos+ leftPos);
    lyCur= topPos - botPos;
    lxCur= rightPos - leftPos;
    
    lxNew=lxCur * exp(-strain); // remember the is a a negative sign in the passed strain for cosmetics only so strain = log(Lcur/Lnew) where conventionaly it is log(Lnew/Lcur)
    lyNew=lyCur * exp(-strain);

    leftPos= xMid-0.5 * lxNew;
    rightPos= xMid+0.5 * lxNew;
    botPos= yMid-0.5 * lyNew;
    topPos= yMid+0.5 * lyNew;
    
    // Apply an affine deformation to all nodal positions keeping the cell center fixed.
    curPosX = curPosX.array() - xMid;
    curPosX *= 1.0/lxCur*lxNew;
    curPosX = curPosX.array() + xMid;
    
    curPosY = curPosY.array() - yMid;
    curPosY *= 1.0/lyCur*lyNew;
    curPosY = curPosY.array() + yMid;
    
}


void Configuration::affine_axial_shearing(const BaseSysData& baseData, const Parameters& pars, double strain){
    
    std::cout << "**** shearing ****   rate " << pars.deformationRate << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
    
    std::cout << "e1  " << e1 <<std::endl;
    std::cout << "phi  " << phi <<std::endl;
    
    yMid=0.5*(topPos+ botPos);
    xMid=0.5*(rightPos+ leftPos);
    lyCur= topPos - botPos;
    lxCur= rightPos - leftPos;
    
    lxNew=lxCur * exp(strain);
    lyNew=lyCur * exp(-strain);
    
    leftPos= xMid-0.5 * lxNew;
    rightPos= xMid+0.5 * lxNew;
    botPos= yMid-0.5 * lyNew;
    topPos= yMid+0.5 * lyNew;
    
    // Apply an affine deformation to all nodal positions keeping the cell center fixed.
    curPosX = curPosX.array() - xMid;
    curPosX *= (1.0/lxCur)*lxNew;
    curPosX = curPosX.array() + xMid;
    
    curPosY = curPosY.array() - yMid;
    curPosY *= (1.0/lyCur)*lyNew;
    curPosY = curPosY.array() + yMid;
    
    
}

void Configuration::special_localized_deformation(const BaseSysData& baseData, const Parameters& pars, const double& gammaX,const double& gammaY,const std::vector<int>& targetNodes){

    yMid=0.5*(curPosY.maxCoeff()+ curPosY.minCoeff());
    xMid=0.5*(curPosX.maxCoeff()+ curPosX.minCoeff());
    lyCur= curPosY.maxCoeff()-curPosY.minCoeff() ;
    lxCur= curPosX.maxCoeff()-curPosX.minCoeff() ;
    
    lxNew=lxCur;
    lyNew=lyCur ;
    
    leftPos= xMid-0.5 * lxNew;
    rightPos= xMid+0.5 * lxNew;
    botPos= yMid-0.5 * lyNew;
    topPos= yMid+0.5 * lyNew;
    

    curPosX=curPosX.array()-xMid;
    curPosX *= exp(gammaX);
    curPosX=curPosX.array()+xMid;
    curPosY=curPosY.array()-yMid;
    curPosY *= exp(gammaY);
    curPosY=curPosY.array()+yMid;
    

}

void Configuration::hold(const BaseSysData& baseData, const Parameters& pars)
{
    
    std::cout << "**** holding **** " << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
    
    std::cout << "e1  " << e1 <<std::endl;
    std::cout << "phi  " << phi <<std::endl;

    yMid=0.5*(topPos+ botPos);
    xMid=0.5*(rightPos+ leftPos);
    lyCur= topPos - botPos;
    lxCur= rightPos - leftPos;
    
    lxNew=lxCur;
    lyNew=lyCur ;
    
    leftPos= xMid-0.5 * lxNew;
    rightPos= xMid+0.5 * lxNew;
    botPos= yMid-0.5 * lyNew;
    topPos= yMid+0.5 * lyNew;
    
    // Apply an affine deformation to all nodal positions keeping the cell center fixed.
    curPosX = curPosX.array() - xMid;
    curPosX *= 1.0/lxCur*lxNew;
    curPosX = curPosX.array() + xMid;
    
    curPosY = curPosY.array() - yMid;
    curPosY *= 1.0/lxCur*lxNew;
    curPosY = curPosY.array() + yMid;
}

void Configuration::update_cells_1(const BaseSysData& baseData, const Parameters& pars)
{
    
    int xCell, yCell, cellId, segment0;
    double x, y;

    

    for (int nodeID=0; nodeID < baseData.numSurfaceNodes; nodeID++)
    {
        x = augmentedCurPosX[baseData.flatSurfaceNodes[nodeID]];
        y = augmentedCurPosY[baseData.flatSurfaceNodes[nodeID]];
        xCell = int( floor( (x-(leftPos-pars.imagesMargin*lxNew))/verletCellSizeX) );
        yCell = int( floor( (y-(botPos-pars.imagesMargin*lyNew))/verletCellSizeY) );

    if ( (xCell < 0) || (yCell < 0) || (xCell >= numXCells)  || (yCell >= numYCells) )
    {
            continue;
    }
        
        cellId = numXCells*yCell+xCell+baseData.numSurfaceNodes;
             
        nodesLinkedList[nodeID] = nodesLinkedList[cellId]; // Note that nodeID is its local ordinal number not the golbal node name assigned in the mesh
        nodesLinkedList[cellId] = nodeID;
        
        if (pars.contactMethod=="nts"){
            segment0 = baseData.nodeToSegments[baseData.flatSurfaceNodes[nodeID]][0];
            
 
            
            segmentsLinkedList_1[segment0] = segmentsLinkedList_1[cellId];
            segmentsLinkedList_1[cellId] = segment0;
        }
        
    }
    
}

void Configuration::update_cells_2(const BaseSysData& baseData, const Parameters& pars)
{
    
    int xCell, yCell, cellId_1, cellId_2, segment0, segment1;
    double x, y;
    
    
    
    for (int nodeID=0; nodeID < baseData.numSurfaceNodes; nodeID++)
    {
        x = augmentedCurPosX[baseData.flatSurfaceNodes[nodeID]];
        y = augmentedCurPosY[baseData.flatSurfaceNodes[nodeID]];
        xCell = int( floor( (x-(leftPos-pars.imagesMargin*lxNew))/verletCellSizeX) );
        yCell = int( floor( (y-(botPos-pars.imagesMargin*lyNew))/verletCellSizeY) );
        
        if ( (xCell < 0) || (yCell < 0) || (xCell >= numXCells)  || (yCell >= numYCells) )
        {
            continue;
        }
        
        segment0 = baseData.nodeToSegments[baseData.flatSurfaceNodes[nodeID]][0];
        segment1 = baseData.nodeToSegments[baseData.flatSurfaceNodes[nodeID]][1];
        
        cellId_1 = numXCells*yCell+xCell+baseData.numSurfaceNodes;
        cellId_2 = numXCells*yCell+xCell;
        
        nodesLinkedList[nodeID] = nodesLinkedList[cellId_1]; // Note that nodeID is its local ordinal number not the golbal node name assigned in the mesh
        nodesLinkedList[cellId_1] = nodeID;
        
        if (segmentsLinkedList_2(segment1,0) == -2) {
            segmentsLinkedList_2(segment1,0) = cellsHeads(cellId_2,0);
            segmentsLinkedList_2(segment1,1) = cellsHeads(cellId_2,1); //inheret the column of the segment associated with this cell
            cellsHeads(cellId_2,0) = segment1;
            cellsHeads(cellId_2,1) = 0;  //column 0 of segmentsLinkedList
            segmentsLinkedList_2(segment1,4) = cellId_2;
        }else if(segmentsLinkedList_2(segment1,2) == -2 && segmentsLinkedList_2(segment1,4) != cellId_2) {
            segmentsLinkedList_2(segment1,2) = cellsHeads(cellId_2,0);
            segmentsLinkedList_2(segment1,3) = cellsHeads(cellId_2,1);//inheret the column of the segment associated with this cell
            cellsHeads(cellId_2,0) = segment1;
            cellsHeads(cellId_2,1) = 2;  //column 2 of segmentsLinkedList
            
            segmentsLinkedList_2(segment1,4) = cellId_2;
        }
        
        if (segmentsLinkedList_2(segment0,0) == -2) {
            segmentsLinkedList_2(segment0,0) = cellsHeads(cellId_2,0);
            segmentsLinkedList_2(segment0,1) = cellsHeads(cellId_2,1); //inheret the column of the segment associated with this cell
            cellsHeads(cellId_2,0) = segment0;
            cellsHeads(cellId_2,1) = 0;  //column 0 of segmentsLinkedList
            
            segmentsLinkedList_2(segment0,4) = cellId_2;
        }else if(segmentsLinkedList_2(segment0,2) == -2 && segmentsLinkedList_2(segment0,4)!= cellId_2 ) {
            segmentsLinkedList_2(segment0,2) = cellsHeads(cellId_2,0);
            segmentsLinkedList_2(segment0,3) = cellsHeads(cellId_2,1);//inheret the column of the segment associated with this cell
            cellsHeads(cellId_2,0) = segment0;
            cellsHeads(cellId_2,1) = 2;  //column 2 of segmentsLinkedList
            
            segmentsLinkedList_2(segment0,4) = cellId_2;
        }
        
    }
    
}

void Configuration::check_force_energy_consistency(const BaseSysData& baseData, const Parameters& pars)
{
    for (int nodeID=0; nodeID < baseData.numOriginalNodes; nodeID++)
    {
        float d = 0.00001;
        compute_forces_pbc(baseData, pars, 0,1,0,0);
        double E1 = totalEnergy;
        curPosX(nodeID) += d;
        compute_forces_pbc(baseData, pars, 0, 1,0,0);
        double E2 = totalEnergy;
        consistencyFactorX(nodeID) = (forceX(nodeID))*d/(E1-E2);
        consistencyErrorFactorX(nodeID) = (forceX(nodeID)*d-(E1-E2))/forceY(nodeID);
        curPosX(nodeID) -= d;
        
        compute_forces_pbc(baseData, pars, 0, 1,0,0);
        E1 = totalEnergy;
        curPosY(nodeID) += d;
        compute_forces_pbc(baseData, pars, 0, 1,0,0);
        E2 = totalEnergy;
        consistencyFactorY(nodeID) = (forceY(nodeID)*d)/(E1-E2);
        consistencyErrorFactorY(nodeID) = (forceY(nodeID)*d-(E1-E2))/forceY(nodeID);
        curPosY(nodeID) -= d;
    }
//    std::cout<<"X"<<std::endl;
//    std::cout<<consistencyFactorX<<std::endl;
//    std::cout<<"Y"<<std::endl;
//    std::cout<<consistencyErrorFactorX<<std::endl;
    
}

void Configuration::fill_augmented_Hessian(){
    
    
        
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        tripletList.reserve(Hixjx.nonZeros()+Hixjy.nonZeros()+Hiyjx.nonZeros()+Hiyjy.nonZeros());
       
        
        //First add the Hixjx values, the indices stay the same as in the augmented H later
        for (int k=0; k<Hixjx.outerSize(); ++k){
            for (Eigen::SparseMatrix<double>::InnerIterator it(Hixjx,k); it; ++it)
            {
                tripletList.push_back(T(int(it.row()),int(it.col()),it.value()));// inner index, here it is equal to it.row()
            }
        }
        //Second add the Hixjy values to right of the fitst block of Hixjx, the row indices shift by #of_nodes
        for (int k=0; k<Hixjy.outerSize(); ++k){
           for (Eigen::SparseMatrix<double>::InnerIterator it(Hixjy,k); it; ++it)
           {
               tripletList.push_back(T(int(it.row()),int(it.col()+Hixjy.cols()),it.value()));// inner index, here it is equal to it.row()
           }
        }
    
        //Now Hiyjx with rows shift
        for (int k=0; k<Hiyjx.outerSize(); ++k){
            for (Eigen::SparseMatrix<double>::InnerIterator it(Hiyjx,k); it; ++it)
            {
                tripletList.push_back(T(int(it.row()+Hiyjx.rows()),int(it.col()),it.value()));// inner index, here it is equal to it.row()
            }
        }
        
    //Now Hiyjx with columns shift
    for (int k=0; k<Hiyjy.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(Hiyjy,k); it; ++it)
        {
            tripletList.push_back(T(int(it.row()+Hiyjy.cols()),int(it.col()+Hiyjy.cols()),it.value()));// inner index, here it is equal to it.row()
        }
    }
        
    Hessian.setFromTriplets(tripletList.begin(), tripletList.end());
    
    
}
