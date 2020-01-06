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
        lyCur =  baseData.lxRef;
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
        //Read the number of nodes found in the first line of the input file
        std::string line,a;
        double b,c,d,x,y;
        std::getline(inFile, line);
        std::istringstream split1(line);
        split1 >> a;
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
    
    consistencyFactorX.resize(baseData.numOriginalNodes);
    consistencyFactorY.resize(baseData.numOriginalNodes);
    consistencyErrorFactorY.resize(baseData.numOriginalNodes);
    consistencyErrorFactorX.resize(baseData.numOriginalNodes);
    
    // If slover is FIRE, initiate zero velocity vectors
    if (pars.solver == "FIRE"){
        velocityX.resize(baseData.numOriginalNodes);
        velocityY.resize(baseData.numOriginalNodes);
        prevVelocityX.resize(baseData.numOriginalNodes);
        prevVelocityY.resize(baseData.numOriginalNodes);
        velocityX.fill(0);
        velocityY.fill(0);
        prevVelocityX.fill(0);
        prevVelocityY.fill(0);
    }
}



void Configuration::update_post_processing_data(const BaseSysData& baseData, const Parameters& pars){
    
    maxR = sqrt((forceX.array()*forceX.array()+forceY.array()*forceY.array()).maxCoeff());
    avgR = sqrt((forceX.array()*forceX.array()+forceY.array()*forceY.array()).mean());
    LX = (curPosX.maxCoeff()-curPosX.minCoeff());
    LY = (curPosY.maxCoeff()-curPosY.minCoeff());
    Fh = 0.5*(-wallForceRight.sum()+wallForceLeft.sum());
    Fv = 0.5*(-wallForceTop.sum()+wallForceBottom.sum());
    P1 = 0.5*(Fv/lxNew+Fh/lyNew);
    P2 = ((CstressYY+CstressYX).dot((refArea.array()*areaRatio.array()).matrix()) + (CstressXX+CstressXY).dot((refArea.array()*areaRatio.array()).matrix()))/(2*LX*LY);
    S1 = 0.5*(Fv/lxNew-Fh/lyNew);
    S2 = ((CstressYY+CstressYX).dot((refArea.array()*areaRatio.array()).matrix())- (CstressXX+CstressXY).dot((refArea.array()*areaRatio.array()).matrix()))/(2*LX*LY);
    ex = log(lxNew/baseData.lxRef);
    ey = log(lyNew/baseData.lyRef);
    A = lxNew * lyNew;
    e0 = 0.5*(ex+ey);
    e1 = (ex-ey);
    phi = pars.Ap / A;
    deltaTotEnergy = totalEnergy - prevTotEnergy;
    prevTotEnergy = totalEnergy;
    L2NormResidual = - (forceX.dot(forceX) + forceY.dot(forceY));
    
    shearVirial = ((CstressYY+CstressYX).dot((refArea.array()*areaRatio.array()).matrix())- (CstressXX+CstressXY).dot((refArea.array()*areaRatio.array()).matrix()));
    pressureVirial = ((CstressYY+CstressYX).dot((refArea.array()*areaRatio.array()).matrix())+ (CstressXX+CstressXY).dot((refArea.array()*areaRatio.array()).matrix()));
}

void Configuration::dump_global_data(const Parameters& pars, char mode, char purpose){
    std::string fname ;
    if (purpose=='i'){
        fname = "/data.txt";
    }else if (purpose=='f'){
        fname = "/final_data.txt";
    }
    if (mode=='w'){
        std::ofstream dataFile;
        dataFile.open (pars.outputFolderName+fname);
        dataFile
        <<  "TotalEnergy"  << std::setw(20)
        <<  "shearVirial"  << std::setw(20)
        <<  "pressureVirial"  << std::setw(20)
        <<  "maxResidualF"  << std::setw(20)
        <<  "pressure1"  << std::setw(20)
        <<  "pressure2"  << std::setw(20)
        <<  "shearStress1"  << std::setw(20)
        <<  "shearStress2"  << std::setw(20)
        <<  "boxArea"  << std::setw(20)
        <<  "phi"  << std::setw(20)
        <<  "e0"  << std::setw(20)
        <<  "e1" << std::endl;
        dataFile.close();
    }else if (mode=='a'){
        std::ofstream dataFile;
        dataFile.open (pars.outputFolderName+fname,  std::ios_base::app);
        dataFile << std::setprecision(9)
        <<  totalEnergy  << std::setw(20)
        <<  shearVirial  << std::setw(20)
        <<  pressureVirial  << std::setw(20)
        <<  maxR  << std::setw(20)
        <<  P1  << std::setw(20)
        <<  P2  << std::setw(20)
        <<  S1  << std::setw(20)
        <<  S2  << std::setw(20)
        <<  A  << std::setw(20)
        <<  phi  << std::setw(20)
        <<  e0  << std::setw(20)
        <<  e1 << std::endl;
        dataFile.close();
        
        
    }
}


void Configuration::dump_per_node(const BaseSysData& baseData, const Parameters& pars, int& timeStep){
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    myfile.open (pars.outputFolderName+"/dataPerNode-"+step+".txt");
    myfile << "numNodes" << "\t" << baseData.numOriginalNodes << std::endl;
    myfile << "kTOverOmega" << "\t" << pars.kTOverOmega << std::endl;
    myfile << "phi" << "\t" << phi << std::endl;
    myfile << "timeStep" << "\t" << timeStep << std::endl;
    myfile << "wallsLRBT" << "\t" << leftPos << "\t" << rightPos << "\t" << botPos << "\t" << topPos << std::endl;
    myfile << "wallForceTop" << "\t" << wallForceTop.sum() << "\t" <<"wallForceBot" << "\t" << wallForceBottom.sum() << "\t" << "wallForceRight" << "\t" << wallForceRight.sum() << "\t" << "wallForceLeft" << "\t" << wallForceLeft.sum() <<  std::endl;
    
    myfile
    << "id"  << std::setw(20)
    <<  "x"  << std::setw(20)
    <<  "y"  << std::setw(20)
    <<  "fx"  << std::setw(20)
    <<  "fy"  << std::setw(20)
    <<  "errfx"  << std::setw(20)
    <<  "errfy"  << std::endl;
    
    for (int i=0; i<curPosX.size();i++ ) {
        myfile
        <<  std::setprecision(9)
        << i << std::setw(20)
        << curPosX[i] << std::setw(20)
        << curPosY[i] << std::setw(20)
        << forceX[i] << std::setw(20)
        << forceY[i] << std::setw(20)
        << consistencyFactorX[i] << std::setw(20)
        << consistencyFactorY[i] << std::endl;
    }
    
    myfile.close();

}

void Configuration::dump_per_node_periodic_images_on(const BaseSysData& baseData, const Parameters& pars, int& timeStep){
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    myfile.open (pars.outputFolderName+"/dataPerNodePeriodicImages-"+step+".txt");
    myfile << "numNodes" << "\t" << baseData.numNodes << std::endl;
    myfile << "kTOverOmega" << "\t" << pars.kTOverOmega << std::endl;
    myfile << "phi" << "\t" << phi << std::endl;
    myfile << "timeStep" << "\t" << timeStep << std::endl;
    myfile << "wallsLRBT" << "\t" << leftPos << "\t" << rightPos << "\t" << botPos << "\t" << topPos << std::endl;
    myfile << "cellSizeX" << "\t" << verletCellSizeX << std::endl;
    myfile << "cellSizeY" << "\t" << verletCellSizeY << std::endl;
    myfile << "effectiveLeftPos" << "\t" << leftPos - pars.imagesMargin*lxNew << std::endl;
    myfile << "effectiveBotPos" << "\t" << botPos - pars.imagesMargin*lyNew << std::endl;
    myfile << "numCellsX" << "\t" << numXCells << std::endl;
    myfile << "numCellsY" << "\t" << numYCells << std::endl;
    myfile << "wallForceTop" << "\t" << wallForceTop.sum() << "\t" <<"wallForceBot" << "\t" << wallForceBottom.sum() << "\t" << "wallForceRight" << "\t" << wallForceRight.sum() << "\t" << "wallForceLeft" << "\t" << wallForceLeft.sum() <<  std::endl;
    
    myfile
    << "id"  << std::setw(20)
    <<  "x"  << std::setw(20)
    <<  "y"  << std::endl;
    
    for (int i=0; i<augmentedCurPosX.size();i++ ) {
        myfile
        <<  std::setprecision(9)
        << i << std::setw(20)
        << augmentedCurPosX[i] << std::setw(20)
        << augmentedCurPosY[i] << std::endl;
    }
    
    myfile.close();
    
}


void Configuration::dump_per_ele(const BaseSysData& baseData, const Parameters& pars, int& timeStep){
    
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    myfile.open (pars.outputFolderName+"/dataPerEle-"+step+".txt");
    myfile << "numElements" << "\t" << baseData.numElements << std::endl;
    myfile << "kTOverOmega" << "\t" << pars.kTOverOmega << std::endl;
    myfile << "internalEnergy" << "\t" <<  std::setprecision(9)<< internalEnergy << std::endl;
    myfile << "wallsEnergy" << "\t" <<  std::setprecision(9)<< wallsEnergy << std::endl;
    myfile << "contactsEnergy" << "\t" <<  std::setprecision(9)<< contactsEnergy << std::endl;
    myfile << "totalEnergy" << "\t" <<  std::setprecision(9)<< totalEnergy << std::endl;
    myfile << "phi" << "\t" << phi << std::endl;
    myfile << "timeStep" << "\t" << timeStep << std::endl;
    myfile << "wallsLRBT" << "\t" << leftPos << "\t" << rightPos << "\t" << botPos << "\t" << topPos << std::endl;
    
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
        << swellingPressurePerEle[i] << std::setw(20)
        << elasticEnergyPerEle[i] << std::setw(20)
        << mixingEnergyPerEle[i] << std::setw(20)
        << internalEnergyPerEle[i] << std::endl;
    }
    
    myfile.close();

}


void Configuration::compress(const BaseSysData& baseData, const Parameters& pars, double strain){
    
    std::cout << "compressing    rate " << pars.deformationRate << "   dt " << pars.dt << "   kWall " << pars.HWallStiffness << "   kSkin " << pars.penaltyStiffness << std::endl;
    
    yMid=0.5*(topPos+ botPos);
    xMid=0.5*(rightPos+ leftPos);
    lyCur= topPos - botPos;
    lxCur= rightPos - leftPos;
    
    lxNew=lxCur * (1 - strain);
    lyNew=lyCur * (1 - strain);

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


void Configuration::shear(const BaseSysData& baseData, const Parameters& pars, double strain){
    

    
    yMid=0.5*(topPos+ botPos);
    xMid=0.5*(rightPos+ leftPos);
    lyCur= topPos - botPos;
    lxCur= rightPos - leftPos;
    
    lxNew=lxCur / (1 - strain);
    lyNew=lyCur * (1 - strain);
    
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
    
    std::cout << "shearing ... " << std::endl;
    std::cout << "e1  " << e1 <<std::endl;
    std::cout << "phi  " << phi <<std::endl;
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
    
    std::cout << "holding ... " << std::endl;
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

void Configuration::update_cells(const BaseSysData& baseData, const Parameters& pars)
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
        
        nodesLinkedList[nodeID] = nodesLinkedList[cellId_1];
        nodesLinkedList[cellId_1] = nodeID;
        
        if (segmentsLinkedList(segment1,0) == -2) {
            segmentsLinkedList(segment1,0) = cellsHeads(cellId_2,0);
            segmentsLinkedList(segment1,1) = cellsHeads(cellId_2,1); //inheret the column of the segment associated with this cell
            cellsHeads(cellId_2,0) = segment1;
            cellsHeads(cellId_2,1) = 0;  //column 0 of segmentsLinkedList
            segmentsLinkedList(segment1,4) = cellId_2;
        }else if(segmentsLinkedList(segment1,2) == -2 && segmentsLinkedList(segment1,4) != cellId_2) {
            segmentsLinkedList(segment1,2) = cellsHeads(cellId_2,0);
            segmentsLinkedList(segment1,3) = cellsHeads(cellId_2,1);//inheret the column of the segment associated with this cell
            cellsHeads(cellId_2,0) = segment1;
            cellsHeads(cellId_2,1) = 2;  //column 2 of segmentsLinkedList
            
            segmentsLinkedList(segment1,4) = cellId_2;
        }
        
        if (segmentsLinkedList(segment0,0) == -2) {
            segmentsLinkedList(segment0,0) = cellsHeads(cellId_2,0);
            segmentsLinkedList(segment0,1) = cellsHeads(cellId_2,1); //inheret the column of the segment associated with this cell
            cellsHeads(cellId_2,0) = segment0;
            cellsHeads(cellId_2,1) = 0;  //column 0 of segmentsLinkedList
            
            segmentsLinkedList(segment0,4) = cellId_2;
        }else if(segmentsLinkedList(segment0,2) == -2 && segmentsLinkedList(segment0,4)!= cellId_2 ) {
            segmentsLinkedList(segment0,2) = cellsHeads(cellId_2,0);
            segmentsLinkedList(segment0,3) = cellsHeads(cellId_2,1);//inheret the column of the segment associated with this cell
            cellsHeads(cellId_2,0) = segment0;
            cellsHeads(cellId_2,1) = 2;  //column 2 of segmentsLinkedList
            
            segmentsLinkedList(segment0,4) = cellId_2;
        }
        
    }
    
}
    
void Configuration::check_force_energy_consistency(const BaseSysData& baseData, const Parameters& pars)
{
    for (int nodeID=0; nodeID < baseData.numOriginalNodes; nodeID++)
    {
        float d = 0.0001;
        compute_forces_PBC(baseData, pars, 0,1);
        double E1 = totalEnergy;
        curPosX(nodeID) += d;
        compute_forces_PBC(baseData, pars, 0, 1);
        double E2 = totalEnergy;
        consistencyFactorX(nodeID) = (forceX(nodeID))*d/(E1-E2);
        consistencyErrorFactorX(nodeID) = (forceX(nodeID)*d-(E1-E2))/forceY(nodeID);
        curPosX(nodeID) -= d;
        
        compute_forces_PBC(baseData, pars, 0, 1);
        E1 = totalEnergy;
        curPosY(nodeID) += d;
        compute_forces_PBC(baseData, pars, 0, 1);
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
