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
#include <chrono>
#include <valarray>
#include "Parameters.hpp"
#include "Configuration.hpp"
#include "BaseSysData.hpp"

void Configuration::NTS_interaction(const int& slaveNodeId, const int& segment,const int& masterMesh,  const BaseSysData& baseData, const Parameters& pars)
{

    int node = baseData.flatSurfaceNodes[slaveNodeId];
    int node0 = baseData.surfaceSegments[segment][0];
    int node1 = baseData.surfaceSegments[segment][1];
    
    double xi = augmentedCurPosX[node];
    double yi = augmentedCurPosY[node];
    double x0 = augmentedCurPosX[node0];
    double y0 = augmentedCurPosY[node0];
    double x1 = augmentedCurPosX[node1];
    double y1 = augmentedCurPosY[node1];

    double dx = x1-x0;
    double dy = y1-y0;
    double L = sqrt(std::pow(dx,2)+std::pow(dy,2));
    double nx = dy/L;
    double ny = - dx/L;
    double signedGap = (xi-x0)*nx+(yi-y0)*ny;
    double s = (xi-x0)*dx/std::pow(L,2)+(yi-y0)*dy/std::pow(L,2);
    double gapSign = 1;
    if (signedGap<0)
    {
       gapSign = -1.0;
    }
    if (s>=0 && s<=1){
        
        if(surNodes_masters(slaveNodeId,0)==-1){
            surNodes_masters(slaveNodeId,0) = masterMesh;
            surNodes_masters(slaveNodeId,1) = segment;
            surNodes_masters(slaveNodeId,2) = 2;
            surNodes_masters(slaveNodeId,3) = signedGap;
        }else if (surNodes_masters(slaveNodeId,0) == masterMesh){
            
            if (abs(surNodes_masters(slaveNodeId,3)) > abs(signedGap))
            {
                surNodes_masters(slaveNodeId,1) = segment;
                surNodes_masters(slaveNodeId,2) = 2;
                surNodes_masters(slaveNodeId,3) = signedGap;
            }
        }else if(surNodes_masters(slaveNodeId,4)==-1){
            surNodes_masters(slaveNodeId,4) = masterMesh;
            surNodes_masters(slaveNodeId,5) = segment;
            surNodes_masters(slaveNodeId,6) = 2;
            surNodes_masters(slaveNodeId,7) = signedGap;
        }else if (surNodes_masters(slaveNodeId,4) == masterMesh){
            
            if (abs(surNodes_masters(slaveNodeId,7)) > abs(signedGap))
            {
                surNodes_masters(slaveNodeId,5) = segment;
                surNodes_masters(slaveNodeId,6) = 2;
                surNodes_masters(slaveNodeId,7) = signedGap;
            }
        }else if(surNodes_masters(slaveNodeId,8)==-1){
            surNodes_masters(slaveNodeId,8) = masterMesh;
            surNodes_masters(slaveNodeId,9) = segment;
            surNodes_masters(slaveNodeId,10) = 2;
            surNodes_masters(slaveNodeId,11) = signedGap;
        }else if (surNodes_masters(slaveNodeId,8) == masterMesh){
            
            if (abs(surNodes_masters(slaveNodeId,11)) > abs(signedGap))
            {
                surNodes_masters(slaveNodeId,9) = segment;
                surNodes_masters(slaveNodeId,10) = 2;
                surNodes_masters(slaveNodeId,11) = signedGap;
            }
        }
        
        
    }else{
        double r0ix = xi-x0 ;
        double r0iy = yi-y0;
        double r1ix = xi-x1;
        double r1iy = yi-y1;
        double signedGap0 = sqrt( std::pow(r0ix,2) + std::pow(r0iy,2) ) * gapSign;
        double signedGap1 = sqrt( std::pow(r1ix,2) + std::pow(r1iy,2) ) * gapSign;
        
        if(surNodes_masters(slaveNodeId,0)==-1){
            
            surNodes_masters(slaveNodeId,0) = masterMesh;
            surNodes_masters(slaveNodeId,1) = segment;
            
            if ( abs(signedGap0) <= abs(signedGap1))
            {
                surNodes_masters(slaveNodeId,2) = 0;
                surNodes_masters(slaveNodeId,3) = signedGap0;
            }else{
                
                surNodes_masters(slaveNodeId,2) = 1;
                surNodes_masters(slaveNodeId,3) = signedGap1;
            }
            
        }else if (surNodes_masters(slaveNodeId,0) == masterMesh){
            
            if (abs(surNodes_masters(slaveNodeId,3)) > abs(signedGap0))
            {
                surNodes_masters(slaveNodeId,1) = segment;
                if ( abs(signedGap0) <= abs(signedGap1))
                {
                    surNodes_masters(slaveNodeId,2) = 0;
                    surNodes_masters(slaveNodeId,3) = signedGap0;
                }else{
                    
                    surNodes_masters(slaveNodeId,2) = 1;
                    surNodes_masters(slaveNodeId,3) = signedGap1;
                }
            }else if (abs(surNodes_masters(slaveNodeId,3)) > abs(signedGap1)){
                surNodes_masters(slaveNodeId,1) = segment;
                surNodes_masters(slaveNodeId,2) = 1;
                surNodes_masters(slaveNodeId,3) = signedGap1;
            }
        }else if(surNodes_masters(slaveNodeId,4)==-1){
            
            surNodes_masters(slaveNodeId,4) = masterMesh;
            surNodes_masters(slaveNodeId,5) = segment;
            
            if ( abs(signedGap0) <= abs(signedGap1))
            {
                surNodes_masters(slaveNodeId,6) = 0;
                surNodes_masters(slaveNodeId,7) = signedGap0;
            }else{
                
                surNodes_masters(slaveNodeId,6) = 1;
                surNodes_masters(slaveNodeId,7) = signedGap1;
            }
            
        }else if (surNodes_masters(slaveNodeId,4) == masterMesh){
            
            if (abs(surNodes_masters(slaveNodeId,7)) > abs(signedGap0))
            {
                surNodes_masters(slaveNodeId,5) = segment;
                if ( abs(signedGap0) <= abs(signedGap1))
                {
                    surNodes_masters(slaveNodeId,6) = 0;
                    surNodes_masters(slaveNodeId,7) = signedGap0;
                }else{
                    
                    surNodes_masters(slaveNodeId,6) = 1;
                    surNodes_masters(slaveNodeId,7) = signedGap1;
                }
            }else if (abs(surNodes_masters(slaveNodeId,7)) > abs(signedGap1)){
                surNodes_masters(slaveNodeId,5) = segment;
                surNodes_masters(slaveNodeId,6) = 1;
                surNodes_masters(slaveNodeId,7) = signedGap1;
            }
        }else if(surNodes_masters(slaveNodeId,8)==-1){
            
            surNodes_masters(slaveNodeId,8) = masterMesh;
            surNodes_masters(slaveNodeId,9) = segment;
            
            if ( abs(signedGap0) <= abs(signedGap1))
            {
                surNodes_masters(slaveNodeId,10) = 0;
                surNodes_masters(slaveNodeId,11) = signedGap0;
            }else{
                
                surNodes_masters(slaveNodeId,10) = 1;
                surNodes_masters(slaveNodeId,11) = signedGap1;
            }
            
        }else if (surNodes_masters(slaveNodeId,8) == masterMesh){
            
            if (abs(surNodes_masters(slaveNodeId,11)) > abs(signedGap0))
            {
                surNodes_masters(slaveNodeId,9) = segment;
                if ( abs(signedGap0) <= abs(signedGap1))
                {
                    surNodes_masters(slaveNodeId,10) = 0;
                    surNodes_masters(slaveNodeId,11) = signedGap0;
                }else{
                    
                    surNodes_masters(slaveNodeId,10) = 1;
                    surNodes_masters(slaveNodeId,11) = signedGap1;
                }
            }else if (abs(surNodes_masters(slaveNodeId,11)) > abs(signedGap1)){
                surNodes_masters(slaveNodeId,9) = segment;
                surNodes_masters(slaveNodeId,10) = 1;
                surNodes_masters(slaveNodeId,11) = signedGap1;
            }
        }

    }

}
