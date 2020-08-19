
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

void Configuration::apply_ntn_repulsions(const BaseSysData &baseData, const Parameters &pars, bool Hessian, const long &timeStep){
   
    for (int cellYid=0; cellYid<numYCells; cellYid++)
       {
           
           for (int cellXid=0; cellXid<numXCells; cellXid++)
           {
               int cellFlatId = numXCells*cellYid + cellXid + baseData.numSurfaceNodes;
               int iNodeId = nodesLinkedList[cellFlatId]; // Note that slaveNodeId is its local ordinal number in surface nodes not the golbal node name assigned in the mesh
        
               
               while (iNodeId>=0)
               {
                   int iMesh = baseData.nodeToMesh[iNodeId];
                   
                   for (auto const& delta: neighborBinDelta)
                   {
                       int  nCellXid = cellXid + delta.first;
                       int  nCellYid = cellYid + delta.second;
                       
                       if ((nCellXid < 0) || (nCellYid < 0) || (nCellXid == numXCells) || (nCellYid == numYCells))
                       {
                           continue;
                       }
                       
                       
                       int nCellFlatId =  numXCells*nCellYid + nCellXid + baseData.numSurfaceNodes;
                       int jNodeId = nodesLinkedList[nCellFlatId];

                       while (jNodeId>=0)
                       {
                           int jMesh = baseData.nodeToMesh[jNodeId];
                           
                           if (iMesh==jMesh)
                           {
                               jNodeId = nodesLinkedList[jNodeId];
                               continue;
                           }
                           
                           // find the global node id from its surfaceNodesList id
                           int iNode = baseData.flatSurfaceNodes[iNodeId];
                           int jNode = baseData.flatSurfaceNodes[jNodeId];
                           
                           if (iNode>=baseData.numOriginalNodes && jNode>=baseData.numOriginalNodes)
                           {
                               jNodeId = nodesLinkedList[jNodeId];
                               continue;
                           }
                        
                           
                           
                           double dxij = augmentedCurPosX[jNode] - augmentedCurPosX[iNode];
                           double dyij = augmentedCurPosY[jNode] - augmentedCurPosY[iNode];
                           double drijSq = pow(dxij,2)+ pow(dyij,2);
                           double drij = sqrt(drijSq);

                           if (pars.ntnRepulsionMethod=="powerlaw"){
                               
                               if (drij > pars.ntnPLRcutoffOverRadius*pars.ntnRadius){
                                   jNodeId = nodesLinkedList[jNodeId];
                                   continue;
                               }
                               double A = - 6*pow(1.0/pars.ntnPLRcutoffOverRadius,10);
                               double forceij= 0.5*pars.ntnPLEnergy/pars.ntnRadius*12*pow((pars.ntnRadius/drij),13)+A/pars.ntnRadius*pow((pars.ntnRadius/drij),3);
                               double forceXij = forceij*dxij/drij;
                               double forceYij = forceij*dyij/drij;
                               
                               if (iNode<baseData.numOriginalNodes){
                                      forceX(iNode) = forceX(iNode)-forceXij;
                                      forceY(iNode) = forceY(iNode)-forceYij;
                                      surfaceForceX(iNode) = surfaceForceX(iNode)-forceXij;
                                      surfaceForceY(iNode) = surfaceForceY(iNode)-forceYij;
                                   
                                   
                                  }
                              if (jNode<baseData.numOriginalNodes){
                                      forceX(jNode) = forceX(jNode)+forceXij;
                                      forceY(jNode) = forceY(jNode)+forceYij;
                                      surfaceForceX(jNode) = surfaceForceX(jNode)+forceXij;
                                      surfaceForceY(jNode) = surfaceForceY(jNode)+forceYij;
                                  }
                                  
                       
                               nodeIinteractions++ ;
                               contactsEnergy+=0.5*pars.ntnPLEnergy*pow((pars.ntnRadius/drij),12)+0.5*A*pow((pars.ntnRadius/drij),2);
                               KWoodXX += - forceij * drij * dxij/drij * dxij/drij ;
                               KWoodYY += - forceij * drij * dyij/drij * dyij/drij ;
                               
                               if (timeStep % pars.dumpEvery == 0 && pars.identifyAndDumbFacets && (pars.ntnRadius/drij) >= 0.4) {
                                     
                                     facets_ntn[std::make_pair(iMesh,jMesh)].push_back(iNode);
                                     facets_ntn[std::make_pair(iMesh,jMesh)].push_back(jNode);
                                     facets_ntn[std::make_pair(iMesh,jMesh)].push_back(drij);
                                     facets_ntn[std::make_pair(iMesh,jMesh)].push_back(forceij);
                                     facets_ntn[std::make_pair(iMesh,jMesh)].push_back(augmentedCurPosX[iNode]);
                                     facets_ntn[std::make_pair(iMesh,jMesh)].push_back(augmentedCurPosX[jNode]);
                                     facets_ntn[std::make_pair(iMesh,jMesh)].push_back(augmentedCurPosY[iNode]);
                                     facets_ntn[std::make_pair(iMesh,jMesh)].push_back(augmentedCurPosY[jNode]);
                                     facets_ntn[std::make_pair(iMesh,jMesh)].push_back(-forceXij);
                                     facets_ntn[std::make_pair(iMesh,jMesh)].push_back(-forceYij);
                                 }
                               
                           }else if (pars.ntnRepulsionMethod=="harmonic") {
                               double d = 2*pars.ntnRadius - drij;
                              
                               if (d<=0){
                                   jNodeId = nodesLinkedList[jNodeId];
                                   continue;
                               }
                               
                               double forceij = pars.ntnHStiffness*d;
                               double forceXij = forceij*dxij/drij;
                               double forceYij = forceij*dyij/drij;
                               
                               if (iNode<=baseData.numOriginalNodes){
                                          forceX(iNode) = forceX(iNode)-forceXij;
                                          forceY(iNode) = forceY(iNode)-forceYij;
                                          surfaceForceX(iNode) = surfaceForceX(iNode)-forceXij;
                                          surfaceForceY(iNode) = surfaceForceY(iNode)-forceYij;
                                      }
                               if (jNode<=baseData.numOriginalNodes){
                                      forceX(jNode) = forceX(jNode)+forceXij;
                                      forceY(jNode) = forceY(jNode)+forceYij;
                                      surfaceForceX(jNode) = surfaceForceX(jNode)+forceXij;
                                      surfaceForceY(jNode) = surfaceForceY(jNode)+forceYij;
                                }
                                      
                               
                                      nodeIinteractions++ ;
                                      contactsEnergy+=0.5*pars.ntnHStiffness*pow(d,2);
                               
                               if (timeStep % pars.dumpEvery == 0 && pars.identifyAndDumbFacets && (pars.ntnRadius/drij) >= 0.2 ) {
                                   
                                   facets_ntn[std::make_pair(iMesh,jMesh)].push_back(iNode);
                                   facets_ntn[std::make_pair(iMesh,jMesh)].push_back(jNode);
                                   facets_ntn[std::make_pair(iMesh,jMesh)].push_back(drij);
                                   facets_ntn[std::make_pair(iMesh,jMesh)].push_back(forceij);
                                   facets_ntn[std::make_pair(iMesh,jMesh)].push_back(augmentedCurPosX[iNode]);
                                   facets_ntn[std::make_pair(iMesh,jMesh)].push_back(augmentedCurPosX[jNode]);
                                   facets_ntn[std::make_pair(iMesh,jMesh)].push_back(augmentedCurPosY[iNode]);
                                   facets_ntn[std::make_pair(iMesh,jMesh)].push_back(augmentedCurPosY[jNode]);
                                   facets_ntn[std::make_pair(iMesh,jMesh)].push_back(-forceXij);
                                   facets_ntn[std::make_pair(iMesh,jMesh)].push_back(-forceYij);
                               }
                           }
                           
                           
                           
        
                           
                           
                           jNodeId = nodesLinkedList[jNodeId];

                       }
                           
                   }
                   
                   iNodeId = nodesLinkedList[iNodeId];
               }
           }
       }
       
    
    
}
