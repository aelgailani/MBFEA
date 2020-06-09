
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

void Configuration::apply_ghost_nodes_to_node_repulsions(const BaseSysData &baseData, const Parameters &pars, bool Hessian, const long &timeStep){
   
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
                           
                           if (iMesh==jMesh )
                           {
                               jNodeId = nodesLinkedList[jNodeId];
                               continue;
                           }
                           
                           int iNode = baseData.flatSurfaceNodes[iNodeId];
                           int jNode = baseData.flatSurfaceNodes[jNodeId];
                           int segment0 =  baseData.nodeToSegments[jNode][0];
                           
                           int node0 = baseData.surfaceSegments[segment0][0];
                           int node1 = baseData.surfaceSegments[segment0][1];
                           
                          
                           double dx01 = augmentedCurPosX[node1] - augmentedCurPosX[node0];
                           double dy01 = augmentedCurPosY[node1] - augmentedCurPosY[node0];
                           
                          
                           for (double s = 0; s < 1.0; s+=1.0/pars.gntn_NGhostNodes){
                               
                               double xs = augmentedCurPosX[node0] + s * dx01;
                               double ys = augmentedCurPosY[node0] + s * dy01;
                               
                               double dxij = xs - augmentedCurPosX[iNode];
                               double dyij = ys - augmentedCurPosY[iNode];
                               double drijSq = pow(dxij,2)+ pow(dyij,2);
                               double drij = sqrt(drijSq);
                               
                               if (drij > pars.ntnPLRcutoff) continue;
                               
                               double forceij=0.5*pars.ntnPLEnergy/pars.ntnRadius*12*pow((pars.ntnRadius/drij),13);
                               double forceXij = forceij*dxij/drij;
                               double forceYij = forceij*dyij/drij;
           
                               if (iNode<=baseData.numOriginalNodes){
                                   forceX(iNode) = forceX(iNode)-forceXij;
                                   forceY(iNode) = forceY(iNode)-forceYij;
                                   surfaceForceX(iNode) = surfaceForceX(iNode)-forceXij;
                                   surfaceForceY(iNode) = surfaceForceY(iNode)-forceYij;
                               }
                               if (node0<=baseData.numOriginalNodes){
                                   forceX(node0) = forceX(node0)+forceXij * (1-s);
                                   forceY(node0) = forceY(node0)+forceYij * (1-s);
                                   surfaceForceX(node0) = surfaceForceX(node0)+forceXij * (1-s);
                                   surfaceForceY(node0) = surfaceForceY(node0)+forceYij * (1-s);
                                   
                                   forceX(node1) = forceX(node1) + forceXij * s;
                                   forceY(node1) = forceY(node1)+  forceYij * s;
                                   surfaceForceX(node1) = surfaceForceX(node1) + forceXij * s;
                                   surfaceForceY(node1) = surfaceForceY(node1) + forceYij * s;
                               }
                                  
                                    
                                 contactsEnergy+=0.5*pars.ntnPLEnergy*pow((pars.ntnRadius/drij),12);
                               
                           }
                           

                           nodeIinteractions++ ;
                           jNodeId = nodesLinkedList[jNodeId];

                       }
                           
                   }
                   
                   iNodeId = nodesLinkedList[iNodeId];
               }
           }
       }
       
    
    
}
