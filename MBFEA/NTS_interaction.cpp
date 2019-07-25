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

void Configuration::NTS_interaction(const int& slaveNodeId, const int& segment,  const BaseSysData& baseData, const Parameters& pars)
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
        if (abs(surNodes_gap[slaveNodeId]) > abs(signedGap))
        {
            surNodes_mSegment[slaveNodeId] = segment;
            surNodes_mSegmentWhichPart[slaveNodeId] = 2;
            surNodes_gap[slaveNodeId] = signedGap ;
//            std::cout << "signedGap    "<< surNodes_gap[slaveNodeId] << std::endl;
        }
    }else{
        double r0ix = xi-x0 ;
        double r0iy = yi-y0;
        double r1ix = xi-x1;
        double r1iy = yi-y1;
        double signedGap0 = sqrt( std::pow(r0ix,2) + std::pow(r0iy,2) ) * gapSign;
        double signedGap1 = sqrt( std::pow(r1ix,2) + std::pow(r1iy,2) ) * gapSign;
        
        if (abs(surNodes_gap[slaveNodeId])> abs(signedGap0))
        {
            surNodes_mSegment[slaveNodeId] = segment;
            surNodes_mSegmentWhichPart[slaveNodeId] = 0;
            surNodes_gap[slaveNodeId] = signedGap0 ;
//            std::cout << "signedGap0    "<< signedGap0 << std::endl;
        }
        
        if (abs(surNodes_gap[slaveNodeId])> abs(signedGap1))
        {
            surNodes_mSegment[slaveNodeId] = segment;
            surNodes_mSegmentWhichPart[slaveNodeId] = 1;
            surNodes_gap[slaveNodeId] = signedGap1 ;
//            std::cout << "signedGap1    "<< signedGap1 << std::endl;
        }

    }

}
