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

void Configuration::NTS_interaction(const int& node, const int& segment,  const BaseSysData& baseData, const Parameters& pars)
{

    
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
    double s = (xi-x0)*dx/std::pow(L,2)+(yi-y0)*dy/std::pow(L,2);
    double gap = (xi-x0)*nx+(yi-y0)*ny;
    double gapSign = 1;
    if (gap<0)
    {
        gapSign = -1;
    }
    
    if (s>=0 && s<=1)
    {
        double f = pars.penaltyStiffness * abs(gap);

        double fx = f * (nx);
        double fy = f * (ny);
        double f0x = -fx * (1-s);
        double f0y = -fy * (1-s);
        double f1x = -fx * (s);
        double f1y = -fy * (s);
        

        if (closestMaster[0]> abs(gap)){
            std::valarray<double> array ={abs(gap),gapSign,fx,fy,double(node0),f0x,f0y,double(node1),f1x,f1y,1};
            closestMaster.swap(array);
            return;
        };
    }else{
        
        double r0ix = xi-x0 ;
        double r0iy = yi-y0;
        double r1ix = xi-x1;
        double r1iy = yi-y1;
        
        double gap0 = sqrt( std::pow(r0ix,2) + std::pow(r0iy,2) );
        double gap1 = sqrt( std::pow(r1ix,2) + std::pow(r1iy,2) );
        
        double f0ix = -pars.penaltyStiffness  * r0ix;
        double f0iy = -pars.penaltyStiffness  * r0iy;
        double f00x = -f0ix;
        double f00y = -f0iy;
        
        double f1ix = -pars.penaltyStiffness  * r1ix;
        double f1iy = -pars.penaltyStiffness  * r1iy;
        double f11x = -f1ix;
        double f11y = -f1iy;
        
       
        if (closestMaster[0]> gap0){
            std::valarray<double> array ={gap0,gapSign,f0ix,f0iy,double(node0),f00x,f00y,double(node1),0,0,0};
            closestMaster.swap(array);
        }
        
        if (closestMaster[0]> gap1){
            std::valarray<double> array ={gap1,gapSign,f1ix,f1iy,double(node0),0,0,double(node1),f11x,f11y,0};
            closestMaster.swap(array);
            
        }
    }
}
