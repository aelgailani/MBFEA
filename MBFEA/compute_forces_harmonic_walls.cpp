
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


void Configuration::compute_forces_harmonic_walls(const BaseSysData& baseData, const Parameters& pars, const long& timeStep, bool surfaceInteractions, bool updatePBC, bool Hessian)
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
    
//    auto finish1 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed1 = finish1 - start1;
//    std::cout << "elapsed time in defGrad:  " << elapsed1.count() << std::endl;
    
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
    

    if (Hessian){
        calculate_the_hessian(pars);
    }

    augmentedCurPosX = curPosX;
    augmentedCurPosY = curPosY;
    
    //This part must come only after the main Hessian because it will be overwritten
    if (surfaceInteractions){
        compute_surface_forces(baseData,pars,Hessian, timeStep);
    }
    
    //apply hamonic repulsion to boudary nodes
    wallForceTop = pars.HWallStiffness*(0.5*((curPosY.array() - topPos).sign()+1))*(topPos-curPosY.array());
    wallForceBottom = pars.HWallStiffness*(0.5*(sign(botPos-curPosY.array())+1))*(botPos-curPosY.array());
//
    wallForceRight = pars.HWallStiffness*(0.5*(sign(curPosX.array()-rightPos)+1))*(rightPos-curPosX.array());
    wallForceLeft = pars.HWallStiffness*(0.5*(sign(leftPos-curPosX.array())+1))*(leftPos-curPosX.array());

    forceX = forceX + wallForceRight + wallForceLeft ;
    forceY = forceY + wallForceTop + wallForceBottom ;
    
    KWoodXX += 1/pars.HWallStiffness*( wallForceRight.array().pow(2).sum() + wallForceLeft.array().pow(2).sum());
    KWoodYY += 1/pars.HWallStiffness*( wallForceTop.array().pow(2).sum() + wallForceBottom.array().pow(2).sum());
//    std::cout << (curPosY.array() - topPos).sign() << std::endl;
    
    wallsEnergy = 0.5/pars.HWallStiffness*( wallForceTop.array().pow(2).sum() + wallForceBottom.array().pow(2).sum() + wallForceRight.array().pow(2).sum() + wallForceLeft.array().pow(2).sum());
    
    
    contactsEnergy += wallsEnergy;
    internalEnergy = internalEnergyPerEle.dot(refArea);
    totalEnergy= internalEnergy + contactsEnergy;
    
//    std::cout << "segmentIinteractions  " << segmentIinteractions << std::endl;
//    std::cout << "nodeIinteractions  " << nodeIinteractions << std::endl;
    std::cout << "interactions  " << nodeIinteractions + segmentIinteractions << std::endl;
//    std::cout << "max skin interference   " << maxInterference << std::endl;
//    std::cout << "max wall interference   " << maxWallinterference << std::endl;
//    std::cout << "min J  " <<areaRatio.array().minCoeff() << std::endl;
    
}

