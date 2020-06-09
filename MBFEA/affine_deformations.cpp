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

void Configuration::affine_axial_shearing_triWalls(const BaseSysData& baseData, const Parameters& pars, double strain, double ctrX, double ctrY){
    
    std::cout << "**** shearing tirangle ****   rate " << pars.deformationRate << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
    std::cout << "strain  " << strain <<std::endl;
//    std::cout << "e1  " << e1 <<std::endl;
//    std::cout << "phi  " << phi <<std::endl;
    
    triAy -=ctrY;
    triAy *= exp(-strain);
    triAy+=ctrY;
    
    triBy -=ctrY;
    triBy *= exp(-strain);
    triBy +=ctrY;
    triCy = triBy;
    
    triBx -=ctrX;
    triBx *= exp(strain);
    triBx +=ctrX;
    
    triCx -=ctrX;
    triCx *= exp(strain);
    triCx +=ctrX;
    
    lyCur= topPos - botPos;
    lxCur= rightPos - leftPos;
    
    topPos = triAy;
    botPos = triBy;
    rightPos = triCx;
    leftPos = triBx;
    
    lyNew= topPos - botPos;
    lxNew= rightPos - leftPos;

    
    // Apply an affine deformation to all nodal positions keeping the cell center fixed.
    curPosX = curPosX.array() - ctrX;
    curPosX *= exp(strain);
    curPosX = curPosX.array() + ctrX;
    
    curPosY = curPosY.array() - ctrY;
    curPosY *= exp(-strain);
    curPosY = curPosY.array() + ctrY;
    
    
}


void Configuration::affine_compression_triWalls(const BaseSysData& baseData, const Parameters& pars, double strain, double ctrX, double ctrY){
    
    std::cout << "**** compressing tirangle ****   rate " << pars.deformationRate << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
    
    std::cout << "strain  " << strain <<std::endl;
//    std::cout << "phi  " << phi <<std::endl;
    
    triAy -=ctrY;
    triAy *= exp(strain);
    triAy+=ctrY;
    
    triBy -=ctrY;
    triBy *= exp(strain);
    triBy +=ctrY;
    triCy = triBy;
    
    triBx -=ctrX;
    triBx *= exp(strain);
    triBx +=ctrX;
    
    triCx -=ctrX;
    triCx *= exp(strain);
    triCx +=ctrX;
    
    lyCur= topPos - botPos;
    lxCur= rightPos - leftPos;
    
    std::cout << "lyCur  " << lyCur <<std::endl;
    std::cout << "lxCur  " << lxCur <<std::endl;
    
    std::cout << "lyCur *= exp(strain) " << lyCur*exp(strain)<<std::endl;
    std::cout << "lxCur *= exp(strain) " << lxCur*exp(strain) <<std::endl;
    
    topPos = triAy;
    botPos = triBy;
    rightPos = triCx;
    leftPos = triBx;
    
    lyNew= topPos - botPos;
    lxNew= rightPos - leftPos;

    std::cout << "lyNew  " << lyNew <<std::endl;
    std::cout << "lxNew  " << lxNew <<std::endl;
    
    
    
    
    // Apply an affine deformation to all nodal positions keeping the cell center fixed.
    curPosX = curPosX.array() - ctrX;
    curPosX *= exp(strain);
    curPosX = curPosX.array() + ctrX;
    
    curPosY = curPosY.array() - ctrY;
    curPosY *= exp(strain);
    curPosY = curPosY.array() + ctrY;
    
    
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
