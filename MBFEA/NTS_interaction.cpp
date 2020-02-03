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
        
        if(surNodes_mMesh1[slaveNodeId]==-1){
            surNodes_mMesh1[slaveNodeId] = masterMesh;
            surNodes_mSegment1[slaveNodeId] = segment;
            surNodes_mPart1[slaveNodeId] = 2;
            surNodes_gap1[slaveNodeId] = signedGap;
            
        }else if (surNodes_mMesh1[slaveNodeId] == masterMesh){
            
            if (abs(surNodes_gap1[slaveNodeId]) > abs(signedGap))
            {
                surNodes_mSegment1[slaveNodeId] = segment;
                surNodes_mPart1[slaveNodeId] = 2;
                surNodes_gap1[slaveNodeId] = signedGap;
            }
        }else if(surNodes_mMesh2[slaveNodeId] ==-1){
            surNodes_mMesh2[slaveNodeId] = masterMesh;
            surNodes_mSegment2[slaveNodeId] = segment;
            surNodes_mPart2[slaveNodeId] = 2;
            surNodes_gap2[slaveNodeId] = signedGap;
            
        }else if (surNodes_mMesh2[slaveNodeId] == masterMesh){
            
            if (abs(surNodes_gap2[slaveNodeId]) > abs(signedGap))
            {
                surNodes_mSegment2[slaveNodeId] = segment;
                surNodes_mPart2[slaveNodeId] = 2;
                surNodes_gap2[slaveNodeId] = signedGap;
            }
        }else if(surNodes_mMesh3[slaveNodeId]==-1){
            surNodes_mMesh3[slaveNodeId] = masterMesh;
            surNodes_mSegment3[slaveNodeId] = segment;
            surNodes_mPart3[slaveNodeId] = 2;
            surNodes_gap3[slaveNodeId] = signedGap;
            
        }else if (surNodes_mMesh3[slaveNodeId] == masterMesh){
            
            if (abs(surNodes_gap3[slaveNodeId]) > abs(signedGap))
            {
                surNodes_mSegment3[slaveNodeId] = segment;
                surNodes_mPart3[slaveNodeId] = 2;
                surNodes_gap3[slaveNodeId] = signedGap;
            }
        }
        
        
    }else{
        double r0ix = xi-x0 ;
        double r0iy = yi-y0;
        double r1ix = xi-x1;
        double r1iy = yi-y1;
        double signedGap0 = sqrt( std::pow(r0ix,2) + std::pow(r0iy,2) ) * gapSign;
        double signedGap1 = sqrt( std::pow(r1ix,2) + std::pow(r1iy,2) ) * gapSign;
        
        if(surNodes_mMesh1[slaveNodeId]==-1){
            
            surNodes_mMesh1[slaveNodeId] = masterMesh;
            surNodes_mSegment1[slaveNodeId] = segment;
            
            if ( abs(signedGap0) <= abs(signedGap1))
            {
                surNodes_mPart1[slaveNodeId] = 0;
                surNodes_gap1[slaveNodeId] = signedGap0;
            }else{
                
                surNodes_mPart1[slaveNodeId] = 1;
                surNodes_gap1[slaveNodeId] = signedGap1;
            }
            
        }else if (surNodes_mMesh1[slaveNodeId] == masterMesh){
            
            if (abs(surNodes_gap1[slaveNodeId]) > abs(signedGap0))
            {
                surNodes_mSegment1[slaveNodeId] = segment;
                if ( abs(signedGap0) <= abs(signedGap1))
                {
                    surNodes_mPart1[slaveNodeId] = 0;
                    surNodes_gap1[slaveNodeId] = signedGap0;
                }else{
                    
                    surNodes_mPart1[slaveNodeId] = 1;
                    surNodes_gap1[slaveNodeId] = signedGap1;
                }
            }else if (abs(surNodes_gap1[slaveNodeId]) > abs(signedGap1)){
                surNodes_mSegment1[slaveNodeId] = segment;
                surNodes_mPart1[slaveNodeId] = 1;
                surNodes_gap1[slaveNodeId] = signedGap1;
            }
        }else if(surNodes_mMesh2[slaveNodeId]==-1){
            
            surNodes_mMesh2[slaveNodeId] = masterMesh;
            surNodes_mSegment2[slaveNodeId] = segment;
            
            if ( abs(signedGap0) <= abs(signedGap1))
            {
                surNodes_mPart2[slaveNodeId] = 0;
                surNodes_gap2[slaveNodeId] = signedGap0;
            }else{
                
                surNodes_mPart2[slaveNodeId] = 1;
                surNodes_gap2[slaveNodeId] = signedGap1;
            }
            
        }else if (surNodes_mMesh2[slaveNodeId] == masterMesh){
            
            if (abs(surNodes_gap2[slaveNodeId]) > abs(signedGap0))
            {
                surNodes_mSegment2[slaveNodeId] = segment;
                if ( abs(signedGap0) <= abs(signedGap1))
                {
                    surNodes_mPart2[slaveNodeId] = 0;
                    surNodes_gap2[slaveNodeId] = signedGap0;
                }else{
                    
                    surNodes_mPart2[slaveNodeId] = 1;
                    surNodes_gap2[slaveNodeId] = signedGap1;
                }
            }else if (abs(surNodes_gap2[slaveNodeId]) > abs(signedGap1)){
                surNodes_mSegment2[slaveNodeId] = segment;
                surNodes_mPart2[slaveNodeId] = 1;
                surNodes_gap2[slaveNodeId] = signedGap1;
            }
        }else if(surNodes_mMesh3[slaveNodeId]==-1){
            
            surNodes_mMesh3[slaveNodeId] = masterMesh;
            surNodes_mSegment3[slaveNodeId] = segment;
            
            if ( abs(signedGap0) <= abs(signedGap1))
            {
                surNodes_mPart3[slaveNodeId] = 0;
                surNodes_gap3[slaveNodeId] = signedGap0;
            }else{
                
                surNodes_mPart3[slaveNodeId] = 1;
                surNodes_gap3[slaveNodeId] = signedGap1;
            }
            
        }else if (surNodes_mMesh3[slaveNodeId] == masterMesh){
            
            if (abs(surNodes_gap3[slaveNodeId]) > abs(signedGap0))
            {
                surNodes_mSegment3[slaveNodeId] = segment;
                if ( abs(signedGap0) <= abs(signedGap1))
                {
                    surNodes_mPart3[slaveNodeId] = 0;
                    surNodes_gap3[slaveNodeId] = signedGap0;
                }else{
                    
                    surNodes_mPart3[slaveNodeId] = 1;
                    surNodes_gap3[slaveNodeId] = signedGap1;
                }
            }else if (abs(surNodes_gap3[slaveNodeId]) > abs(signedGap1)){
                surNodes_mSegment3[slaveNodeId] = segment;
                surNodes_mPart3[slaveNodeId] = 1;
                surNodes_gap3[slaveNodeId] = signedGap1;
            }
        }

    }

}
