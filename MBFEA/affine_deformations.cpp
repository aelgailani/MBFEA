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


void Configuration::affine_compression(const BaseSysData& baseData, const Parameters& pars, double strain, double ctrX, double ctrY, long timestep){
    
    if(timestep%pars.writeToConsoleEvery==0){
        std::cout << "**** compressing ****   rate " << pars.deformationRate << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
    }
    
    lyCur= topPos - botPos;
    lxCur= rightPos - leftPos;
    
    lxNew=lxCur * exp(-strain); // remember the is a a negative sign in the passed strain for cosmetics only so strain = log(Lcur/Lnew) where conventionaly it is log(Lnew/Lcur)
    lyNew=lyCur * exp(-strain);

    leftPos= ctrX-0.5 * lxNew;
    rightPos= ctrX+0.5 * lxNew;
    botPos= ctrY-0.5 * lyNew;
    topPos= ctrY+0.5 * lyNew;
    
    // Apply an affine deformation to all nodal positions keeping the cell center fixed.
      // Apply an affine deformation to all nodal positions keeping the cell center fixed.
      curPosX = curPosX.array() - ctrX;
      if(timestep%pars.dumpEvery) homVx = (curPosX*(exp(strain)-1))/ pars.dt;
      curPosX *= exp(strain);
      curPosX = curPosX.array() + ctrX;
      
      curPosY = curPosY.array() - ctrY;
      if(timestep%pars.dumpEvery) homVy = (curPosY*(exp(strain)-1))/ pars.dt;
      curPosY *= exp(strain);
      curPosY = curPosY.array() + ctrY;
    
}


void Configuration::affine_axial_shearing(const BaseSysData& baseData, const Parameters& pars, double strain, double ctrX, double ctrY, long timestep){
    
    if(timestep%pars.writeToConsoleEvery==0){
        std::cout << "**** shearing ****   rate " << pars.deformationRate << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
        
        std::cout << "e1  " << e1 <<std::endl;
        std::cout << "phi  " << phi <<std::endl;
    }

    lyCur= topPos - botPos;
    lxCur= rightPos - leftPos;
    
    lxNew=lxCur * exp(strain);
    lyNew=lyCur * exp(-strain);
    
    leftPos= ctrX-0.5 * lxNew;
    rightPos= ctrX+0.5 * lxNew;
    botPos= ctrY-0.5 * lyNew;
    topPos= ctrY+0.5 * lyNew;
    
    // Apply an affine deformation to all nodal positions keeping the cell center fixed.
    // Apply an affine deformation to all nodal positions keeping the cell center fixed.
    curPosX = curPosX.array() - ctrX;
    if(timestep%pars.dumpEvery) homVx = (curPosX*(exp(strain)-1))/ pars.dt;
    curPosX *= exp(strain);
    curPosX = curPosX.array() + ctrX;
    
    curPosY = curPosY.array() - ctrY;
    if(timestep%pars.dumpEvery) homVy = (curPosY*(exp(-strain)-1))/ pars.dt;
    curPosY *= exp(-strain);
    curPosY = curPosY.array() + ctrY;
    
    
}

void Configuration::affine_axial_shearing_triWalls(const BaseSysData& baseData, const Parameters& pars, double strain, double ctrX, double ctrY, long timestep){
    
     if(timestep%pars.writeToConsoleEvery==0){
        std::cout << "**** shearing tirangle ****   rate " << pars.deformationRate << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
        std::cout << "strain  " << strain <<std::endl;
    //    std::cout << "e1  " << e1 <<std::endl;
    //    std::cout << "phi  " << phi <<std::endl;
     }
    height*=exp(-strain);
    base*=exp(strain);
    
    triAy = ctrY+height*2.0/3.0;
    triBy = ctrY-height*1.0/3.0;
    triCy = triBy;
    
    triBx =ctrX-0.5*base;
    triCx =ctrX+0.5*base;

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
    if(timestep%pars.dumpEvery) homVx = (curPosX*(exp(strain)-1))/ pars.dt;
    curPosX *= exp(strain);
    curPosX = curPosX.array() + ctrX;
    
    curPosY = curPosY.array() - ctrY;
    if(timestep%pars.dumpEvery) homVy = (curPosY*(exp(-strain)-1))/ pars.dt;
    curPosY *= exp(-strain);
    curPosY = curPosY.array() + ctrY;
    
    
}


void Configuration::affine_compression_triWalls(const BaseSysData& baseData, const Parameters& pars, double strain, double ctrX, double ctrY, long timestep){
    
     if(timestep%pars.writeToConsoleEvery==0){
        std::cout << "**** compressing tirangle ****   rate " << pars.deformationRate << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
        
        std::cout << "strain  " << strain <<std::endl;
    //    std::cout << "phi  " << phi <<std::endl;
     }
    height*=exp(strain);
    base*=exp(strain);
    
    triAy = ctrY+height*2.0/3.0;
    triBy = ctrY-height*1.0/3.0;
    triCy = triBy;
    
    triBx =ctrX-0.5*base;
    triCx =ctrX+0.5*base;
    
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
    if(timestep%pars.dumpEvery) homVx = (curPosX*(exp(strain)-1))/ pars.dt;
    curPosX *= exp(strain);
    curPosX = curPosX.array() + ctrX;
    
    curPosY = curPosY.array() - ctrY;
    if(timestep%pars.dumpEvery) homVy = (curPosY*(exp(strain)-1))/ pars.dt;
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
