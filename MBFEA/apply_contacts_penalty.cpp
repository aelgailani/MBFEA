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


void Configuration::apply_contacts_penalty(const BaseSysData& baseData, const Parameters& pars, const std::valarray<int>& surNodes_mMesh, const std::valarray<int>& surNodes_mSegment, const std::valarray<int>& surNodes_mPart, const std::valarray<double>& surNodes_gap, bool Hessian)
{

    for (int nodeID = 0; nodeID <baseData.numSurfaceNodes; nodeID++){
            
                if(surNodes_mMesh[nodeID]==-1){
                    continue;
                }
                
                int whichPart =  surNodes_mPart[nodeID];
                double gap = surNodes_gap[nodeID];
                // first check if there is interference
                if (gap>= 0)
                {
                    continue;
                }
                if (abs(gap) > maxInterference)
                {
                    maxInterference = abs(gap);
                }
            
                //now the rest comes
            
                int node = baseData.flatSurfaceNodes[nodeID];
                int segment = surNodes_mSegment[nodeID];
                double xi = augmentedCurPosX[node];
                double yi = augmentedCurPosY[node];
                int node0 = baseData.surfaceSegments[segment][0];
                int node1 = baseData.surfaceSegments[segment][1];
            
            
                if (whichPart==2){
                    

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
                    double f = pars.penaltyStiffness * abs(gap);
                    double fx = f * (nx);
                    double fy = f * (ny);
                    double f0x = -fx * (1-s);
                    double f0y = -fy * (1-s);
                    double f1x = -fx * (s);
                    double f1y = -fy * (s);

                    if (node < baseData.numOriginalNodes)
                    {
                        forceX(node) = forceX(node) + fx ;
                        forceY(node) = forceY(node) + fy ;
                        surfaceForceX(node) = surfaceForceX(node) + fx; // the surfaceForces vectors are just for debugging purposes here
                        surfaceForceY(node) = surfaceForceY(node) + fy;
                        
                    }
                    if ( node0 < baseData.numOriginalNodes){
                        forceX(node0) = forceX(node0) + f0x ;
                        forceY(node0) = forceY(node0) + f0y ;
                        surfaceForceX(node0) = surfaceForceX(node0) + f0x;
                        surfaceForceY(node0) = surfaceForceY(node0) + f0y;
                    }
                    if ( node1 < baseData.numOriginalNodes){
                        forceX(node1) = forceX(node1) + f1x ;
                        forceY(node1) = forceY(node1) + f1y ;
                        surfaceForceX(node1) = surfaceForceX(node1) + f1x;
                        surfaceForceY(node1) = surfaceForceY(node1) + f1y;
                        
                    }
                    
                    if ( node < baseData.numOriginalNodes || node0 < baseData.numOriginalNodes){
                    contactsEnergy += pars.penaltyStiffness/2 *(gap*gap);
                    segmentIinteractions++ ;
                    }
                    
                    // add the Hessian part
                    if (Hessian) {
                        add_d1_contributions_to_Hessian(pars.penaltyStiffness, xi,yi,x0,y0,x1,y1, node, node0, node1, baseData);
                    }
                    
                    
                }else if(whichPart==0){
                    
                    double x0 = augmentedCurPosX[node0];
                    double y0 = augmentedCurPosY[node0];
                    double r0ix = xi-x0;
                    double r0iy = yi-y0;
                    double f0ix = -pars.penaltyStiffness  * r0ix;
                    double f0iy = -pars.penaltyStiffness  * r0iy;
                    double f00x = -f0ix;
                    double f00y = -f0iy;
                    
                    if (node < baseData.numOriginalNodes)
                    {
                        forceX(node) = forceX(node) + f0ix ;
                        forceY(node) = forceY(node) + f0iy ;
                        surfaceForceX(node) = surfaceForceX(node) + f0ix;
                        surfaceForceY(node) = surfaceForceY(node) + f0iy;
                    }
                    if ( node0 < baseData.numOriginalNodes){
                        
                        forceX(node0) = forceX(node0) + f00x ;
                        forceY(node0) = forceY(node0) + f00y ;
                        surfaceForceX(node0) = surfaceForceX(node0) + f00x;
                        surfaceForceY(node0) = surfaceForceY(node0) + f00y;
                    }
                    
                    if ( node < baseData.numOriginalNodes || node0 < baseData.numOriginalNodes){
                        contactsEnergy += pars.penaltyStiffness/2 *(gap*gap);
                        segmentIinteractions++ ;
                        
                    }
                     
                    if (Hessian) {
                    add_d2_contributions_to_Hessian(pars.penaltyStiffness, xi,yi,x0,y0, node, node0, baseData);
                     }
//
                }else if(whichPart==1){
                    
                    double x1 = augmentedCurPosX[node1];
                    double y1 = augmentedCurPosY[node1];
                    double r1ix = xi-x1;
                    double r1iy = yi-y1;
                    
                    double f1ix = -pars.penaltyStiffness  * r1ix;
                    double f1iy = -pars.penaltyStiffness  * r1iy;
                    double f11x = -f1ix;
                    double f11y = -f1iy;
                    
                    if (node < baseData.numOriginalNodes)
                    {
                        forceX(node) = forceX(node) + f1ix ;
                        forceY(node) = forceY(node) + f1iy ;
                        surfaceForceX(node) = surfaceForceX(node) + f1ix;
                        surfaceForceY(node) = surfaceForceY(node) + f1iy;
                    }
                    if ( node1 < baseData.numOriginalNodes){
                        
                        forceX(node1) = forceX(node1) + f11x ;
                        forceY(node1) = forceY(node1) + f11y ;
                        surfaceForceX(node1) = surfaceForceX(node1) + f11x;
                        surfaceForceY(node1) = surfaceForceY(node1) + f11y;
                    }
                    if ( node < baseData.numOriginalNodes || node1 < baseData.numOriginalNodes){
                        contactsEnergy += pars.penaltyStiffness/2 *(gap*gap);
                        segmentIinteractions++ ;
                        
                    }
                    
                    if (Hessian){
                    add_d2_contributions_to_Hessian(pars.penaltyStiffness, xi,yi,x1,y1, node, node1, baseData);
                    }
                }
            
        }
}
