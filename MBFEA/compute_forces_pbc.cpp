
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
//#include <chrono>
#include <valarray>
#include "Parameters.hpp"
#include "Configuration.hpp"
#include "BaseSysData.hpp"


void Configuration::compute_forces_pbc(const BaseSysData& baseData, const Parameters& pars, const long& timeStep, bool surfaceInteractions, bool updatePBC, bool Hessian)
{
//    auto start1 = std::chrono::high_resolution_clock::now();
    
    //erase previous step data
    contactsEnergy=0;
    KWoodYY=0;
    KWoodXX=0;
    // initialize surface interaction forces on nodes to zero
    surfaceForceX.resize(baseData.numOriginalNodes);
    surfaceForceY.resize(baseData.numOriginalNodes);
    surfaceForceX.fill(0);
    surfaceForceY.fill(0);
    //compute the deformation gradient
    defGradXX = gradX * curPosX;
    defGradYX = gradX * curPosY;
    defGradXY = gradY * curPosX;
    defGradYY = gradY * curPosY;
    
    
    //compute the determinant of defGrad
    areaRatio = defGradXX.array() * defGradYY.array() - defGradXY.array() * defGradYX.array();

    
//    auto finish2 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed2 = finish2 - finish1;
//    std::cout << "elapsed time in areaRatio:  " << elapsed2.count() << std::endl;
    
    //compute the magnitude squared of defGrad
    fSquared = defGradXX.array().pow(2)+defGradXY.array().pow(2)+defGradYX.array().pow(2)+defGradYY.array().pow(2);
    
//    auto finish3 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed3 = finish3 - finish2;
//    std::cout << "elapsed time in fSquared:  " << elapsed3.count() << std::endl;
    
    //compute the local elastic energy density
    elasticEnergyPerEle = pars.NkT/2.0 * (fSquared.array()- 2 - 2*log(abs(areaRatio.array())));
//    auto finish4 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed4 = finish4 - finish3;
//    std::cout << "elapsed time in elasticEnergyPerEle:  " << elapsed4.count() << std::endl;
    
    //compute the local mixing energy density
    mixingEnergyPerEle = pars.kTOverOmega*(areaRatio.array()-1.0)*(log((abs(areaRatio.array())-1.0)/abs(areaRatio.array()))+pars.chi/abs(areaRatio.array()));
//    auto finish5 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed5 = finish5 - finish4;
//    std::cout << "elapsed time in mixingEnergyPerEle:  " << elapsed5.count() << std::endl;
    
    //compute the local total energy density
    internalEnergyPerEle = elasticEnergyPerEle.array() + mixingEnergyPerEle.array();
//    auto finish6 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed6 = finish6 - finish5;
//    std::cout << "elapsed time in internalEnergyPerEle:  " << elapsed6.count() << std::endl;
    
    //compute the inverse of defGrad
    invDefGradTransXX = defGradYY.array()/ areaRatio.array();
    invDefGradTransXY = defGradXY.array()/ areaRatio.array() * -1 ;
    invDefGradTransYX = defGradYX.array()/ areaRatio.array() * -1;
    invDefGradTransYY = defGradXX.array()/ areaRatio.array();
    
//    auto finish7 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed7 = finish7 - finish6;
//    std::cout << "elapsed time in invDefGradTrans:  " << elapsed7.count() << std::endl;
    
    swellingPressurePerEle = pars.kTOverOmega * ( (pars.chi+areaRatio.array()) / (areaRatio.array()).pow(2) + log((areaRatio.array()-1)/areaRatio.array()) );
    
    PK1stressXX = pars.NkT * defGradXX.array() + ( swellingPressurePerEle.array() * areaRatio.array() - pars.NkT ) * invDefGradTransXX.array();
    PK1stressXY = pars.NkT * defGradXY.array() + ( swellingPressurePerEle.array() * areaRatio.array() - pars.NkT ) * invDefGradTransYX.array();
    PK1stressYX = pars.NkT * defGradYX.array() + ( swellingPressurePerEle.array() * areaRatio.array() - pars.NkT ) * invDefGradTransXY.array();
    PK1stressYY = pars.NkT * defGradYY.array() + ( swellingPressurePerEle.array() * areaRatio.array() - pars.NkT ) * invDefGradTransYY.array();
    
    CstressXX = pars.NkT/areaRatio.array()*(defGradXX.array()*defGradXX.array()+defGradXY.array()*defGradXY.array())+(swellingPressurePerEle.array()-1/areaRatio.array())*(defGradXX.array()*invDefGradTransXX.array()+defGradXY.array()*invDefGradTransYX.array());
    CstressXY = pars.NkT/areaRatio.array()*(defGradXX.array()*defGradYX.array()+defGradXY.array()*defGradYY.array())+(swellingPressurePerEle.array()-1/areaRatio.array())*(defGradYX.array()*invDefGradTransXX.array()+defGradYY.array()*invDefGradTransYX.array());
    CstressYX = pars.NkT/areaRatio.array()*(defGradYX.array()*defGradXX.array()+defGradYY.array()*defGradXY.array())+(swellingPressurePerEle.array()-1/areaRatio.array())*(defGradXX.array()*invDefGradTransXY.array()+defGradXY.array()*invDefGradTransYY.array());
    CstressYY = pars.NkT/areaRatio.array()*(defGradYX.array()*defGradYX.array()+defGradYY.array()*defGradYY.array())+(swellingPressurePerEle.array()-1/areaRatio.array())*(defGradYX.array()*invDefGradTransXY.array()+defGradYY.array()*invDefGradTransYY.array());
    
//    auto finish11 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed11 = finish11 - finish7;
//    std::cout << "elapsed time in stresses:  " << elapsed11.count() << std::endl;
//
    
    
    //compute the nodal forces from the stresses
    forceX = - gradX.transpose() * (PK1stressXX.array() * abs(refArea.array())).matrix() - gradY.transpose() * (PK1stressXY.array() * abs(refArea.array())).matrix();
    forceY = - gradX.transpose() * (PK1stressYX.array() * abs(refArea.array())).matrix() - gradY.transpose() * (PK1stressYY.array() * abs(refArea.array())).matrix();
    interForceX = forceX;
    interForceY = forceY;
    
//    auto finish12 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed12 = finish12 - finish11;
//    std::cout << "elapsed time in internal forces:  " << elapsed12.count() << std::endl;
    
    // Create images x y
    if (updatePBC) {
        curPosXL = curPosX.array()-lxNew;
        curPosXR = curPosX.array()+lxNew;
        curPosXB = curPosX;
        curPosXT = curPosX;
        curPosXBL = curPosX.array()-lxNew;
        curPosXBR = curPosX.array()+lxNew;
        curPosXTL = curPosX.array()-lxNew;
        curPosXTR = curPosX.array()+lxNew;
        
        curPosYL = curPosY;
        curPosYR = curPosY;
        curPosYB = curPosY.array()-lyNew;
        curPosYT = curPosY.array()+lyNew;
        curPosYBL = curPosY.array()-lyNew;
        curPosYBR = curPosY.array()-lyNew;
        curPosYTL = curPosY.array()+lyNew;
        curPosYTR = curPosY.array()+lyNew;

    }
    
    augmentedCurPosX.resize(baseData.numNodes);
    augmentedCurPosY.resize(baseData.numNodes);
    augmentedCurPosX << curPosX, curPosXL, curPosXR, curPosXB, curPosXT, curPosXBL, curPosXBR, curPosXTL, curPosXTR;
    augmentedCurPosY << curPosY, curPosYL, curPosYR, curPosYB, curPosYT, curPosYBL, curPosYBR, curPosYTL, curPosYTR;
//    auto finish13 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed13 = finish13 - finish12;
//    std::cout << "elapsed time in periodic images:  " << elapsed13.count() << std::endl;
    
    if (Hessian){
        calculate_the_hessian(pars);
    }
//    auto finish14 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed14 = finish14 - finish13;
//    std::cout << "elapsed time in augmentedCurPos:  " << elapsed14.count() << std::endl;

//    auto finish15 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed15 = finish15 - start1;
//    std::cout << "Total time elapsed in internal nodes:  " << elapsed15.count() << std::endl;
    
    //This part must come only after the main Hessian because it will be overwritten
    if (surfaceInteractions){
        compute_surface_forces(baseData,pars,Hessian, timeStep);
    }
    
//    auto finish16 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed16 = finish16 - finish15;
//    std::cout << "Total time elapsed in surface interactions:  " << elapsed16.count() << std::endl;
    
    internalEnergy = internalEnergyPerEle.dot(refArea);
    totalEnergy= internalEnergy + contactsEnergy;
    
//    std::cout << "segmentIinteractions  " << segmentIinteractions << std::endl;
//    std::cout << "nodeIinteractions  " << nodeIinteractions << std::endl;
    std::cout << "interactions  " << nodeIinteractions + segmentIinteractions << std::endl;
//    std::cout << "max skin interference   " << maxInterference << std::endl;
//    std::cout << "max wall interference   " << maxWallinterference << std::endl;
//    std::cout << "min J  " <<areaRatio.array().minCoeff() << std::endl;
    
    

    
//    Prints for debugging purposes
    
//    std::cout << "gradX.transpose()  \n" <<gradX.transpose()<< std::endl;
//    std::cout << "gradY.transpose()  \n" <<gradY.transpose()<< std::endl;
//    std::cout << "Fxx \n " <<defGradXX<< std::endl;
//    std::cout << "Fxy^\n " <<defGradXY<< std::endl;
//    std::cout << "Fyx \n " <<defGradYX<< std::endl;
//    std::cout << "Fyy \n " <<defGradYY<< std::endl;
//    std::cout << "Fxx^-1 \n " <<invDefGradTransXX<< std::endl;
//    std::cout << "Fxy^-1 \n " <<invDefGradTransXY<< std::endl;
//    std::cout << "Fyx^-1 \n " <<invDefGradTransYX<< std::endl;
//    std::cout << "Fyy^-1 \n " <<invDefGradTransYY<< std::endl;
//    std::cout << "refAreas \n " <<refArea<< std::endl;
//    std::cout << "J \n " <<areaRatio<< std::endl;

    
//    std::cout << "W''J^2 \n " <<WmPrimePrimeJSquared<< std::endl;
//    std::cout << "W'J \n " <<swellingPressurePerEle.array()*areaRatio.array()<< std::endl;
//    std::cout << "Kxxxx \n " <<Kxxxx<< std::endl;
//    std::cout << "Kxxxy \n " <<Kxxxy<< std::endl;
//    std::cout << "Kxxyx \n " <<Kxxyx<< std::endl;
//    std::cout << "Kxxyy \n " <<Kxxyy<< std::endl;
//    std::cout << "Kxyxy \n " <<Kxyxy<< std::endl;
//    std::cout << "Kxyyx \n " <<Kxyyx<< std::endl;
//    std::cout << "Kxyyy \n " <<Kxyyy<< std::endl;
//    std::cout << "Kyxyx \n " <<Kyxyx<< std::endl;
//    std::cout << "Kyxyy \n " <<Kyxyy<< std::endl;
//    std::cout << "Kyyyy \n " <<Kyyyy<< std::endl;
//
//
//    std::cout << "KMxxjx \n " <<KMxxjx<< std::endl;
//    std::cout << "KMxxjy \n " <<KMxxjy<< std::endl;
//    std::cout << "KMxyjx \n " <<KMxyjx<< std::endl;
//    std::cout << "KMxyjy \n " <<KMxyjy<< std::endl;
//    std::cout << "KMyxjx \n " <<KMyxjx<< std::endl;
//    std::cout << "KMyxjy \n " <<KMyxjy<< std::endl;
//    std::cout << "KMyyjx \n " <<KMyyjx<< std::endl;
//    std::cout << "KMyyjy \n " <<KMyyjy<< std::endl;
    
//    std::cout << "Hixjx \n " <<Hixjx<< std::endl;
//    std::cout << "Hixjy \n " <<Hixjy<< std::endl;
//    std::cout << "Hiyjx \n " <<Hiyjx<< std::endl;
//    std::cout << "Hiyjy \n " <<Hiyjy<< std::endl;
    
//    std::cout << "gradX.col(0) \n" <<gradX.col(0).cwiseProduct(Kxxxx) << std::endl;
//    std::cout << "gradX.col(1) \n" <<gradX.col(1) << std::endl;
//    std::cout << "Hixjx " <<Hixjx << std::endl;
//    std::cout << "Hixjy " <<Hixjx << std::endl;

}
