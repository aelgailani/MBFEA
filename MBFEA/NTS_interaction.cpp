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

    int masterMesh = baseData.surfaceSegments[segment][2];
    
    int node0 = baseData.surfaceSegments[segment][0];
    int node1 = baseData.surfaceSegments[segment][1];
    
    double xi = augmentedCurPosX[node];
    double yi = augmentedCurPosY[node];
    double x0 = augmentedCurPosX[node0];
    double y0 = augmentedCurPosY[node0];
    double x1 = augmentedCurPosX[node1];
    double y1 = augmentedCurPosY[node1];
//    double xi = curPosX[node];
//    double yi = curPosY[node];
//    double x0 = curPosX[node0];
//    double y0 = curPosY[node0];
//    double x1 = curPosX[node1];
//    double y1 = curPosY[node1];
// 
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
        double f0 = -f*(1-s);
        double f1 = -f * s;
        double fx = f * (nx);
        double fy = f * (ny);
        double f0x = -fx * (1-s);
        double f0y = -fy * (1-s);
        double f1x = -fx * (s);
        double f1y = -fy * (s);
        
        if (gaps[{node,masterMesh}].size()==0) {
            gaps[{node,masterMesh}]={abs(gap),gapSign,double(node),f,fx,fy,double(node0),f0,f0x,f0y,double(node1),f1,f1x,f1y,nx,ny,s,double(segment),(xi-x0)*nx,(yi-y0)*ny };
            return;
        }
        if (gaps[{node,masterMesh}][0]>abs(gap)){
            gaps[{node,masterMesh}] = {abs(gap),gapSign,double(node),f,fx,fy,double(node0),f0,f0x,f0y,double(node1),f1,f1x,f1y,nx,ny,s,double(segment),(xi-x0)*nx,(yi-y0)*ny };
            return;
        };
    }else{
        
        double r0ix = xi-x0 ;
        double r0iy = yi-y0;
        double r1ix = xi-x1;
        double r1iy = yi-y1;
        
        double gap0 = sqrt( std::pow(r0ix,2) + std::pow(r0iy,2) );
        double gap1 = sqrt( std::pow(r1ix,2) + std::pow(r1iy,2) );
        
        double f0i= pars.penaltyStiffness  * gap0;
        double f0ix = -pars.penaltyStiffness  * r0ix;
        double f0iy = -pars.penaltyStiffness  * r0iy;
        double f00= -f0i;
        double f00x = -f0ix;
        double f00y = -f0iy;
        
        double f1i = pars.penaltyStiffness  * gap1;
        double f1ix = -pars.penaltyStiffness  * r1ix;
        double f1iy = -pars.penaltyStiffness  * r1iy;
        double f11 = -f1i;
        double f11x = -f1ix;
        double f11y = -f1iy;
        
        if (gaps[{node,masterMesh}].size()==0){
            gaps[{node,masterMesh}] = {gap0,gapSign,double(node),f0i,f0ix,f0iy,double(node0),f00,f00x,f00y,double(segment),r0ix,r0iy};
            
        }
        
        if (gaps[{node,masterMesh}][0]>gap0){
            gaps[{node,masterMesh}] = {gap0,gapSign,double(node),f0i,f0ix,f0iy,double(node0),f00,f00x,f00y,double(segment),r0ix,r0iy};
            
        }
        
        if (gaps[{node,masterMesh}][0]>gap1){
            gaps[{node,masterMesh}] = {gap1,gapSign,double(node),f1i,f1ix,f1iy,double(node1),f11,f11x,f11y,double(segment),r1ix,r1iy};
            
        }
    }
}
