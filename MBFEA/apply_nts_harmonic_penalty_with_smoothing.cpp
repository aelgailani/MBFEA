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
//#include <chrono
#include <valarray>
#include "Parameters.hpp"
#include "Configuration.hpp"
#include "BaseSysData.hpp"
#include "utility_functions.hpp"

void Configuration::apply_nts_harmonic_penalty_with_smoothing(const BaseSysData& baseData, const Parameters& pars, const std::valarray<int>& surNodes_mMesh, const std::valarray<int>& surNodes_mSegment, const std::valarray<int>& surNodes_mPart, const std::valarray<double>& surNodes_gap, const std::valarray<double>& surNodes_smoothCurveX, const std::valarray<double>& surNodes_smoothCurveY,const std::valarray<double>& surNodes_smoothCurve_g, bool Hessian, const long& timeStep)
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
            
               double x0 = augmentedCurPosX[node0];
               double y0 = augmentedCurPosY[node0];
               double x1 = augmentedCurPosX[node1];
               double y1 = augmentedCurPosY[node1];
                
                if (whichPart==2){
                    

                   
                    double dx = x1-x0;
                    double dy = y1-y0;
                    double L = sqrt(std::pow(dx,2)+std::pow(dy,2));
                    double nx = dy/L;
                    double ny = - dx/L;
                    double s = (xi-x0)*dx/std::pow(L,2)+(yi-y0)*dy/std::pow(L,2);
                    double f = pars.ntsHarmonicPenaltyStiffness * abs(gap);
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
                        contactsEnergy += pars.ntsHarmonicPenaltyStiffness/2 *(gap*gap);
                        segmentIinteractions++ ;
                            
                        double Dd1Dxi=(y0 - y1)/pow(pow(-x0 + x1,2) + pow(-y0 + y1,2),0.5);

                        double Dd1Dyi=(-x0 + x1)/pow(pow(-x0 + x1,2) + pow(-y0 + y1,2),0.5);

                        double Dd1Dx0=(1.*(x0 - 1.*x1)*(x0 - 1.*xi)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (-y0 + y1)/pow(pow(x0 - x1,2) + pow(y0 - y1,2),0.5) + (y0 - yi)/pow(pow(x0 - x1,2) + pow(y0 - y1,2),0.5) - (1.*pow(x0 - 1.*x1,2)*(y0 - 1.*yi))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

                        double Dd1Dy0=(x0 - x1)/pow(pow(x0 - x1,2) + pow(y0 - y1,2),0.5) + (-x0 + xi)/pow(pow(x0 - x1,2) + pow(y0 - y1,2),0.5) + (1.*(x0 - 1.*xi)*pow(y0 - 1.*y1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) - (1.*(x0 - 1.*x1)*(y0 - 1.*y1)*(y0 - 1.*yi))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

                        double Dd1Dx1=(-1.*(x0 - 1.*x1)*(x0 - 1.*xi)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (1.*pow(x0 - 1.*x1,2)*(y0 - 1.*yi))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (-y0 + yi)/pow(pow(x0 - x1,2) + pow(y0 - y1,2),0.5);

                        double Dd1Dy1=(x0 - xi)/pow(pow(x0 - x1,2) + pow(y0 - y1,2),0.5) - (1.*(x0 - 1.*xi)*pow(y0 - 1.*y1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (1.*(x0 - 1.*x1)*(y0 - 1.*y1)*(y0 - 1.*yi))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

                            KWoodXX += -f*(Dd1Dxi*xi+Dd1Dx0*x0+Dd1Dx1*x1);
                            KWoodYY += -f*(Dd1Dyi*yi+Dd1Dy0*y0+Dd1Dy1*y1);
//                        std::cout << KWoodXX << std::endl;
//                        std::cout << KWoodYY << std::endl;
                        if(timeStep % pars.dumpEvery == 0 && pars.dumpSmoothenCurves){
                            int masterMesh = surNodes_mMesh[nodeID];
                            interactions_nts[std::make_pair(node,masterMesh)].push_back(xi);
                            interactions_nts[std::make_pair(node,masterMesh)].push_back(yi);
                            interactions_nts[std::make_pair(node,masterMesh)].push_back(x0+s*dx);
                            interactions_nts[std::make_pair(node,masterMesh)].push_back(y0+s*dy);
                        }
                    }
                    
                    // add the Hessian part
                    if (Hessian) {
                        add_d1_contributions_to_the_hessian(pars.ntsHarmonicPenaltyStiffness, xi,yi,x0,y0,x1,y1, node, node0, node1, baseData);
                    }
                    
                    //add the facets elements if required
                    if (timeStep % pars.dumpEvery == 0 && pars.identifyAndDumbFacets) {
                        int masterMesh = surNodes_mMesh[nodeID];
                        int slaveMesh = baseData.nodeToSegments[node][2];
                        facets[std::make_pair(slaveMesh,masterMesh)].push_back(node);
                        facets[std::make_pair(slaveMesh,masterMesh)].push_back(node0);
                        facets[std::make_pair(slaveMesh,masterMesh)].push_back(node1);
                    }
                    
                    
                    
                }else if(whichPart==1){
                    int node2 =  baseData.surfaceSegments[baseData.surfaceSegments[segment][3]][0];
                    double x2= augmentedCurPosX[node2];
                    double y2= augmentedCurPosY[node2];
                    double xg = surNodes_smoothCurveX[nodeID];
                    double yg = surNodes_smoothCurveY[nodeID];
                    double g = surNodes_smoothCurve_g[nodeID];
                    double Nnorm = sqrt((xg-xi)*(xg-xi)+(yg-yi)*(yg-yi));
                    double nx = (xg - xi)/Nnorm;  ///// because g always point in negative direction of the normal here
                    double ny = (yg - yi)/Nnorm;
                    
                    double DxDx2 = pars.alpha_HermitPol/4.0*(g*g-2*g+1);
                    double DyDy2 = DxDx2;
                    double DxDx0 = 1/4.*(4-2*pars.alpha_HermitPol*(g*g+1));
                    double DyDy0 = DxDx0;
                    double DxDx1 = abs(gap) * pars.alpha_HermitPol/4.0*(g*g+2*g+1);
                    double DyDy1 = DxDx1;
        
                    double f1 = pars.ntsHarmonicPenaltyStiffness * abs(gap);
                    double f2 = pars.ntsHarmonicPenaltyStiffness * abs(gap) * DyDy2;
                    double f3 = pars.ntsHarmonicPenaltyStiffness * abs(gap)*DyDy0;
                    double f4 = pars.ntsHarmonicPenaltyStiffness * abs(gap) * DyDy1;

   
                    
                    if (node < baseData.numOriginalNodes)
                    {
                        forceX(node) = forceX(node) + f1*nx;
                        forceY(node) = forceY(node) + f1*ny;
                        surfaceForceX(node) = surfaceForceX(node) + f1*nx;
                        surfaceForceY(node) = surfaceForceY(node) + f1*ny;
                    }
                    
                    if ( node0 < baseData.numOriginalNodes){
                        
                        forceX(node0) = forceX(node0) - f3*nx ;
                        forceY(node0) = forceY(node0) - f3*ny ;
                        surfaceForceX(node0) = surfaceForceX(node0) - f3*nx;
                        surfaceForceY(node0) = surfaceForceY(node0) - f3*ny;
                    }
                    
                    if ( node1 < baseData.numOriginalNodes){
                        
                        forceX(node1) = forceX(node1) - f4*nx ;
                        forceY(node1) = forceY(node1) - f4*ny ;
                        surfaceForceX(node1) = surfaceForceX(node1) - f4*nx;
                        surfaceForceY(node1) = surfaceForceY(node1) - f4*ny;
                    }
                    
                    if ( node2 < baseData.numOriginalNodes){
                        
                        forceX(node2) = forceX(node2) - f2*nx ;
                        forceY(node2) = forceY(node2) - f2*ny ;
                        surfaceForceX(node2) = surfaceForceX(node2) - f2*nx;
                        surfaceForceY(node2) = surfaceForceY(node2) - f2*ny;
                    }
                    
                    if ( node < baseData.numOriginalNodes || node0 < baseData.numOriginalNodes || node2 < baseData.numOriginalNodes){
                        contactsEnergy += pars.ntsHarmonicPenaltyStiffness/2 *(gap*gap);
                        nodeIinteractions++ ;
                        
                            ////add Kirkwood piece
                        double kd = pars.ntsHarmonicPenaltyStiffness*abs(gap);
                        
                        
                        double DdDxi =(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*x0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*x1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*x2)/4. + xi)/sqrt(pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*x0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*x1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*x2)/4. + xi,2) + pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*y0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*y1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*y2)/4. + yi,2));
                        double DdDyi =(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*y0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*y1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*y2)/4. + yi)/sqrt(pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*x0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*x1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*x2)/4. + xi,2) + pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*y0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*y1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*y2)/4. + yi,2));
                        double DdDx0 =((-4 + 2*pars.alpha_HermitPol*(1 + pow(g,2)))*(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*x0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*x1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*x2)/4. + xi))/(4.*sqrt(pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*x0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*x1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*x2)/4. + xi,2) + pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*y0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*y1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*y2)/4. + yi,2)));
                        double DdDy0 =((-4 + 2*pars.alpha_HermitPol*(1 + pow(g,2)))*(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*y0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*y1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*y2)/4. + yi))/(4.*sqrt(pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*x0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*x1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*x2)/4. + xi,2) + pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*y0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*y1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*y2)/4. + yi,2)));
                        double DdDx1 =-(pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*x0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*x1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*x2)/4. + xi))/(4.*sqrt(pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*x0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*x1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*x2)/4. + xi,2) + pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*y0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*y1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*y2)/4. + yi,2)));
                        double DdDy1 =-(pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*y0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*y1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*y2)/4. + yi))/(4.*sqrt(pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*x0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*x1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*x2)/4. + xi,2) + pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*y0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*y1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*y2)/4. + yi,2)));
                        double DdDx2 =-(pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*x0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*x1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*x2)/4. + xi))/(4.*sqrt(pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*x0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*x1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*x2)/4. + xi,2) + pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*y0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*y1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*y2)/4. + yi,2)));
                        double DdDy2 =-(pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*y0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*y1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*y2)/4. + yi))/(4.*sqrt(pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*x0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*x1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*x2)/4. + xi,2) + pow(-((4 - 2*pars.alpha_HermitPol*(1 + pow(g,2)))*y0)/4. - (pars.alpha_HermitPol*(1 + 2*g + pow(g,2))*y1)/4. - (pars.alpha_HermitPol*(1 - 2*g + pow(g,2))*y2)/4. + yi,2)));
                        
                        KWoodXX +=-kd*(DdDxi*xi+DdDx0*x0+DdDx1*x1+DdDx2*x2);
                        KWoodXX +=-kd*(DdDyi*yi+DdDy0*y0+DdDy1*y1+DdDy2*y2);
                        

                        if(timeStep % pars.dumpEvery == 0 && pars.dumpSmoothenCurves){
                            int masterMesh = surNodes_mMesh[nodeID];
                            interactions_nts[std::make_pair(node,masterMesh)].push_back(xi);
                            interactions_nts[std::make_pair(node,masterMesh)].push_back(yi);
                            interactions_nts[std::make_pair(node,masterMesh)].push_back(xg);
                            interactions_nts[std::make_pair(node,masterMesh)].push_back(yg);
                        }
                    }
                    
                    
                    
//
//                    if (Hessian) {
//                    add_d2_contributions_to_the_hessian(pars.ntsHarmonicPenaltyStiffness, xi,yi,x0,y0, node, node0, baseData);
//                     }
                    
                    //add the facets elements if required
                    if (timeStep % pars.dumpEvery == 0 && pars.identifyAndDumbFacets) {
                        int masterMesh = surNodes_mMesh[nodeID];
                        int slaveMesh = baseData.nodeToSegments[node][2];
                        facets[std::make_pair(slaveMesh,masterMesh)].push_back(node);
                        facets[std::make_pair(slaveMesh,masterMesh)].push_back(node0);
                    }
//
                }
            
        }
}

