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
    
    curLY_W= topPos - botPos;
    curLX_W= rightPos - leftPos;
    
    newLX_W=curLX_W * exp(-strain); // remember the is a a negative sign in the passed strain for cosmetics only so strain = log(Lcur/Lnew) where conventionaly it is log(Lnew/Lcur)
    newLY_W=curLY_W * exp(-strain);

    leftPos= ctrX-0.5 * newLX_W;
    rightPos= ctrX+0.5 * newLX_W;
    botPos= ctrY-0.5 * newLY_W;
    topPos= ctrY+0.5 * newLY_W;
    
    // Apply an affine deformation to all nodal positions keeping the cell center fixed.
      // Apply an affine deformation to all nodal positions keeping the cell center fixed.
      curPosX = curPosX.array() - ctrX;
      if(timestep%pars.dumpEvery) homVx = (curPosX*(exp(-strain)-1))/ pars.dt;
      curPosX *= exp(-strain);
      curPosX = curPosX.array() + ctrX;
      
      curPosY = curPosY.array() - ctrY;
      if(timestep%pars.dumpEvery) homVy = (curPosY*(exp(-strain)-1))/ pars.dt;
      curPosY *= exp(-strain);
      curPosY = curPosY.array() + ctrY;
    
}


void Configuration::affine_axial_shearing(const BaseSysData& baseData, const Parameters& pars, double strain, double ctrX, double ctrY, long timestep){
    
    if(timestep%pars.writeToConsoleEvery==0){
        std::cout << "**** shearing ****   rate " << pars.deformationRate << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
        
        std::cout << "e1  " << e1_W <<std::endl;
        std::cout << "phi  " << phi <<std::endl;
    }

    curLY_W= topPos - botPos;
    curLX_W= rightPos - leftPos;
    
    newLX_W=curLX_W * exp(strain);
    newLY_W=curLY_W * exp(-strain);
    
    leftPos= ctrX-0.5 * newLX_W;
    rightPos= ctrX+0.5 * newLX_W;
    botPos= ctrY-0.5 * newLY_W;
    topPos= ctrY+0.5 * newLY_W;
    
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

void Configuration::surface_nodes_affine_axial_shearing(const BaseSysData& baseData, const Parameters& pars, double strain, double ctrX, double ctrY, long timestep){
    
    if(timestep%pars.writeToConsoleEvery==0){
        std::cout << "**** shearing surface nodes only ****   rate " << pars.deformationRate << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
        
        std::cout << "e1  " << e1_W <<std::endl;
        std::cout << "phi  " << phi <<std::endl;
    }

    curLY_W= topPos - botPos;
    curLX_W= rightPos - leftPos;
    
    newLX_W=curLX_W * exp(strain);
    newLY_W=curLY_W * exp(-strain);
    
    leftPos= ctrX-0.5 * newLX_W;
    rightPos= ctrX+0.5 * newLX_W;
    botPos= ctrY-0.5 * newLY_W;
    topPos= ctrY+0.5 * newLY_W;
    
    // Apply an affine deformation to surface nodal positions keeping the cell center fixed.
    // Apply an affine deformation to surface nodal positions keeping the cell center fixed.
    homVx = curPosX*0;
    homVy = curPosX*0;
    for (int nodeID=0; nodeID < baseData.numOriginalSurfaceNodes; nodeID++)
    {
           curPosX[baseData.flatSurfaceNodes[nodeID]] = curPosX[baseData.flatSurfaceNodes[nodeID]] - ctrX;
           if(timestep%pars.dumpEvery) homVx[baseData.flatSurfaceNodes[nodeID]] = (curPosX[baseData.flatSurfaceNodes[nodeID]]*(exp(strain)-1))/ pars.dt;
           curPosX[baseData.flatSurfaceNodes[nodeID]] = curPosX[baseData.flatSurfaceNodes[nodeID]] * exp(strain);
           curPosX[baseData.flatSurfaceNodes[nodeID]] = curPosX[baseData.flatSurfaceNodes[nodeID]] + ctrX;
           
           curPosY[baseData.flatSurfaceNodes[nodeID]] =curPosY[baseData.flatSurfaceNodes[nodeID]] - ctrY;
           if(timestep%pars.dumpEvery) homVy[baseData.flatSurfaceNodes[nodeID]] = (curPosY[baseData.flatSurfaceNodes[nodeID]]*(exp(-strain)-1))/ pars.dt;
           curPosY[baseData.flatSurfaceNodes[nodeID]] =curPosY[baseData.flatSurfaceNodes[nodeID]]* exp(-strain);
           curPosY[baseData.flatSurfaceNodes[nodeID]] = curPosY[baseData.flatSurfaceNodes[nodeID]] + ctrY;
           
        
    }
          
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

    curLY_W= topPos - botPos;
    curLX_W= rightPos - leftPos;
    
    topPos = triAy;
    botPos = triBy;
    rightPos = triCx;
    leftPos = triBx;
    
    newLY_W= topPos - botPos;
    newLX_W= rightPos - leftPos;

    
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
    
    curLY_W= topPos - botPos;
    curLX_W= rightPos - leftPos;

    
    topPos = triAy;
    botPos = triBy;
    rightPos = triCx;
    leftPos = triBx;
    
    newLY_W= topPos - botPos;
    newLX_W= rightPos - leftPos;

    
    
    
    
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
    curLY_W= curPosY.maxCoeff()-curPosY.minCoeff() ;
    curLX_W= curPosX.maxCoeff()-curPosX.minCoeff() ;
    
    newLX_W=curLX_W;
    newLY_W=curLY_W ;
    
    leftPos= xMid-0.5 * newLX_W;
    rightPos= xMid+0.5 * newLX_W;
    botPos= yMid-0.5 * newLY_W;
    topPos= yMid+0.5 * newLY_W;
    

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
    
    std::cout << "e1  " << e1_W <<std::endl;
    std::cout << "phi  " << phi <<std::endl;

    yMid=0.5*(topPos+ botPos);
    xMid=0.5*(rightPos+ leftPos);
    curLY_W= topPos - botPos;
    curLX_W= rightPos - leftPos;
    
    newLX_W=curLX_W;
    newLY_W=curLY_W ;
    
    leftPos= xMid-0.5 * newLX_W;
    rightPos= xMid+0.5 * newLX_W;
    botPos= yMid-0.5 * newLY_W;
    topPos= yMid+0.5 * newLY_W;
    
    // Apply an affine deformation to all nodal positions keeping the cell center fixed.
    curPosX = curPosX.array() - xMid;
    curPosX *= 1.0/curLX_W*newLX_W;
    curPosX = curPosX.array() + xMid;
    
    curPosY = curPosY.array() - yMid;
    curPosY *= 1.0/curLX_W*newLX_W;
    curPosY = curPosY.array() + yMid;
}
