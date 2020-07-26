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
#include "utility_functions.hpp"

void Configuration::nts_find_closest_approach_with_smoothing(const int& slaveNodeId, const int& segment, const int& masterMesh,  const BaseSysData& baseData, const Parameters& pars)
{

    int node = baseData.flatSurfaceNodes[slaveNodeId];
    int node0 = baseData.surfaceSegments[segment][0];
    int node1 = baseData.surfaceSegments[segment][1];
    
//    std::cout<< "x1  " << node << std::endl;
//    std::cout<< "x2  " << node0 << std::endl;
//    std::cout<< "x3  " << node1 << std::endl;
    
    
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
    
    
    
    
    if (s>=pars.alpha_HermitPol && s<=(1-pars.alpha_HermitPol)){
        
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
        
        
    }
//    else if (s>=(1-pars.alpha_HermitPol)){///// in case s lies right to the segment then check the downstream curve. Notice here x2-->node0, x3-->node1, x4-->node2
//        int node2 =  baseData.surfaceSegments[baseData.nodeToSegments[node1][1]][1];
//        double x2 = augmentedCurPosX[node2];
//        double y2 = augmentedCurPosY[node2];
//
//        ///// claculate s for the right segemet
//        dx = x2-x1;
//        dy = y2-y1;
//        L = sqrt(std::pow(dx,2)+std::pow(dy,2));
//        nx = dy/L;
//        ny = - dx/L;
//        signedGap = (xi-x1)*nx+(yi-y1)*ny;
//        s = (xi-x1)*dx/std::pow(L,2)+(yi-y1)*dy/std::pow(L,2);
////        std::cout<< "x4  " << node2 << std::endl;
//        if (s>pars.alpha_HermitPol) return; //// This needs to be revised. weather useing return or give a value to the gap and then return
//        struct closestPointOnHermitianCurve point = find_closest_point_on_Hermitian_interpolation(xi, x0, x1, x2, yi, y0, y1, y2, pars.alpha_HermitPol);
//
//
//
//        if(surNodes_mMesh1[slaveNodeId]==-1){
//
//            surNodes_mMesh1[slaveNodeId] = masterMesh;
//            surNodes_mSegment1[slaveNodeId] = segment; // segment , curves have smae ID which is left node (node0) when you walk on the screen leaving the center of the mesh to your left
//            surNodes_mPart1[slaveNodeId] = 1;
//            surNodes_gap1[slaveNodeId] = point.gap;
//            surNodes_smoothCurveX1[slaveNodeId] = point.x;
//            surNodes_smoothCurveY1[slaveNodeId] = point.y;
//            surNodes_smoothCurve_g1[slaveNodeId] = point.g;
//
//        }else if (surNodes_mMesh1[slaveNodeId] == masterMesh){
//
//            if (abs(surNodes_gap1[slaveNodeId]) > abs(point.gap))
//            {
//                surNodes_mSegment1[slaveNodeId] = segment;
//                surNodes_mPart1[slaveNodeId] = 1;
//                surNodes_gap1[slaveNodeId] = point.gap;
//                surNodes_smoothCurveX1[slaveNodeId] = point.x;
//                surNodes_smoothCurveY1[slaveNodeId] = point.y;
//                surNodes_smoothCurve_g1[slaveNodeId] = point.g;
//
//            }
//        }else if(surNodes_mMesh2[slaveNodeId]==-1){
//
//            surNodes_mMesh2[slaveNodeId] = masterMesh;
//            surNodes_mSegment2[slaveNodeId] = segment; // segment , curves have smae ID which is left node (node0) when you walk on the screen leaving the center of the mesh to your left
//
//
//            surNodes_mPart2[slaveNodeId] = 1;
//            surNodes_gap2[slaveNodeId] = point.gap;
//            surNodes_smoothCurveX2[slaveNodeId] = point.x;
//            surNodes_smoothCurveY2[slaveNodeId] = point.y;
//            surNodes_smoothCurve_g2[slaveNodeId] = point.g;
//
//        }else if (surNodes_mMesh2[slaveNodeId] == masterMesh){
//
//            if (abs(surNodes_gap2[slaveNodeId]) > abs(point.gap))
//            {
//                surNodes_mSegment2[slaveNodeId] = segment;
//                surNodes_mPart2[slaveNodeId] = 1;
//                surNodes_gap2[slaveNodeId] = point.gap;
//                surNodes_smoothCurveX2[slaveNodeId] = point.x;
//                surNodes_smoothCurveY2[slaveNodeId] = point.y;
//                surNodes_smoothCurve_g2[slaveNodeId] = point.g;
//
//            }
//        }else if(surNodes_mMesh3[slaveNodeId]==-1){
//
//
//            surNodes_mMesh3[slaveNodeId] = masterMesh;
//            surNodes_mSegment3[slaveNodeId] = segment; // segment , curves have smae ID which is left node (node0) when you walk on the screen leaving the center of the mesh to your left
//
//
//            surNodes_mPart3[slaveNodeId] = 1;
//            surNodes_gap3[slaveNodeId] = point.gap;
//            surNodes_smoothCurveX3[slaveNodeId] = point.x;
//            surNodes_smoothCurveY3[slaveNodeId] = point.y;
//            surNodes_smoothCurve_g3[slaveNodeId] = point.g;
//
//        }else if (surNodes_mMesh3[slaveNodeId] == masterMesh){
//            if (abs(surNodes_gap3[slaveNodeId]) > abs(point.gap))
//            {
//                surNodes_mSegment3[slaveNodeId] = segment;
//                surNodes_mPart3[slaveNodeId] = 1;
//                surNodes_gap3[slaveNodeId] = point.gap;
//                surNodes_smoothCurveX3[slaveNodeId] = point.x;
//                surNodes_smoothCurveY3[slaveNodeId] = point.y;
//                surNodes_smoothCurve_g3[slaveNodeId] = point.g;
//
//            }
//
//        }
//
//    }
    else if (s<pars.alpha_HermitPol){///// in case s lies left to the segment then check the upstream curve. Notice here x2-->node2, x3-->node0, x4-->node1
        
        int node2 =  baseData.surfaceSegments[baseData.surfaceSegments[segment][3]][0]; //// baseData.surfaceSegments[segment][3] gives downstream segment ID
        double x2 = augmentedCurPosX[node2];
        double y2 = augmentedCurPosY[node2];
        
        ///// claculate s for the left segemet
        dx = x0-x2;
        dy = y0-y2;
        L = sqrt(std::pow(dx,2)+std::pow(dy,2));
        nx = dy/L;
        ny = - dx/L;
        signedGap = (xi-x2)*nx+(yi-y2)*ny;
        double s0 = (xi-x2)*dx/std::pow(L,2)+(yi-y2)*dy/std::pow(L,2);
//        std::cout<< "x4  " << node2 << std::endl;
        if (s0<(1-pars.alpha_HermitPol)) return; //// This needs to be revised. weather useing return or give a value to the gap and then return
        
        
        struct closestPointOnHermitianCurve point = find_closest_point_on_Hermitian_interpolation(xi, x2, x0, x1, yi, y2, y0, y1, pars.alpha_HermitPol);
        
     
        if(surNodes_mMesh1[slaveNodeId]==-1){
            
            surNodes_mMesh1[slaveNodeId] = masterMesh;
            surNodes_mSegment1[slaveNodeId] = baseData.surfaceSegments[segment][3]; ///// segment , curves have smae ID which is left node (node0) when you walk on the screen leaving the center of the mesh to your left.  baseData.surfaceSegments[segment][4] gives downstream segment ID
            surNodes_mPart1[slaveNodeId] = 1;
            surNodes_gap1[slaveNodeId] = point.gap;
            surNodes_smoothCurveX1[slaveNodeId] = point.x;
            surNodes_smoothCurveY1[slaveNodeId] = point.y;
            surNodes_smoothCurve_g1[slaveNodeId] = point.g;
            
        }else if (surNodes_mMesh1[slaveNodeId] == masterMesh){
            
            if (abs(surNodes_gap1[slaveNodeId]) > abs(point.gap))
            {
                surNodes_mSegment1[slaveNodeId] = baseData.surfaceSegments[segment][3];
                surNodes_mPart1[slaveNodeId] = 1;
                surNodes_gap1[slaveNodeId] = point.gap;
                surNodes_smoothCurveX1[slaveNodeId] = point.x;
                surNodes_smoothCurveY1[slaveNodeId] = point.y;
                surNodes_smoothCurve_g1[slaveNodeId] = point.g;
                
            }
        }else if(surNodes_mMesh2[slaveNodeId]==-1){
            
            surNodes_mMesh2[slaveNodeId] = masterMesh;
            surNodes_mSegment2[slaveNodeId] = baseData.surfaceSegments[segment][4]; /////////// segment , curves have smae ID which is left node (node0) when you walk on the screen leaving the center of the mesh to your left .   baseData.surfaceSegments[segment][4] gives downstream segment ID
                
            
            surNodes_mPart2[slaveNodeId] = 1;
            surNodes_gap2[slaveNodeId] = point.gap;
            surNodes_smoothCurveX2[slaveNodeId] = point.x;
            surNodes_smoothCurveY2[slaveNodeId] = point.y;
            surNodes_smoothCurve_g2[slaveNodeId] = point.g;
            
        }else if (surNodes_mMesh2[slaveNodeId] == masterMesh){
            
            if (abs(surNodes_gap2[slaveNodeId]) > abs(point.gap))
            {
                surNodes_mSegment2[slaveNodeId] = baseData.surfaceSegments[segment][3];
                surNodes_mPart2[slaveNodeId] = 1;
                surNodes_gap2[slaveNodeId] = point.gap;
                surNodes_smoothCurveX2[slaveNodeId] = point.x;
                surNodes_smoothCurveY2[slaveNodeId] = point.y;
                surNodes_smoothCurve_g2[slaveNodeId] = point.g;
                
            }
        }else if(surNodes_mMesh3[slaveNodeId]==-1){
            
      
            surNodes_mMesh3[slaveNodeId] = masterMesh;
            surNodes_mSegment3[slaveNodeId] = baseData.surfaceSegments[segment][3]; // segment , curves have smae ID which is left node (node0) when you walk on the screen leaving the center of the mesh to your left
                
            
            surNodes_mPart3[slaveNodeId] = 1;
            surNodes_gap3[slaveNodeId] = point.gap;
            surNodes_smoothCurveX3[slaveNodeId] = point.x;
            surNodes_smoothCurveY3[slaveNodeId] = point.y;
            surNodes_smoothCurve_g3[slaveNodeId] = point.g;
            
        }else if (surNodes_mMesh3[slaveNodeId] == masterMesh){
            if (abs(surNodes_gap3[slaveNodeId]) > abs(point.gap))
            {
                surNodes_mSegment3[slaveNodeId] = baseData.surfaceSegments[segment][3];
                surNodes_mPart3[slaveNodeId] = 1;
                surNodes_gap3[slaveNodeId] = point.gap;
                surNodes_smoothCurveX3[slaveNodeId] = point.x;
                surNodes_smoothCurveY3[slaveNodeId] = point.y;
                surNodes_smoothCurve_g3[slaveNodeId] = point.g;
                
            }
            
        }

    }

}
