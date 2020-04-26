
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
                           
                           if (iMesh==jMesh )
                           {
                               jNodeId = nodesLinkedList[jNodeId];
                               continue;
                           }
                           
                           // find the global node id from its surfaceNodesList id
                           int iNode = baseData.flatSurfaceNodes[iNodeId];
                           int jNode = baseData.flatSurfaceNodes[jNodeId];
                           double dxij = augmentedCurPosX[jNode] - augmentedCurPosX[iNode];
                           double dyij = augmentedCurPosY[jNode] - augmentedCurPosY[iNode];
                           double drijSq = pow(dxij,2)+ pow(dyij,2);
                           double drij = sqrt(drijSq);
                           double forceij=0.5*pars.ntnRepulseEnergy/pars.ntnLjScale*12*pow((pars.ntnLjScale/drij),13);
                           
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
                           contactsEnergy+=0.5*pars.ntnRepulseEnergy*pow((pars.ntnLjScale/drij),12);
                           
                           jNodeId = nodesLinkedList[jNodeId];

                       }
                           
                   }
                   
                   iNodeId = nodesLinkedList[iNodeId];
               }
           }
       }
       
    
    
}
