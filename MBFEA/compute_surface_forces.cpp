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

void Configuration::compute_surface_forces(const BaseSysData& baseData, const Parameters& pars, const int& timeStep)
{
    maxInterference = 0;
    segmentIinteractions = 0;
    nodeIinteractions = 0;

    
    numXCells = floor(lxNew*(1+2*pars.imagesMargin)/pars.verletCellCutoff);
    numYCells = floor(lyNew*(1+2*pars.imagesMargin)/pars.verletCellCutoff);
    verletCellSizeX  = lxNew*(1+2*pars.imagesMargin)/numXCells;
    verletCellSizeY = lyNew*(1+2*pars.imagesMargin)/numYCells;
    nodesLinkedList.resize(numXCells*numYCells+baseData.numSurfaceNodes,-1);
//    surNodes_gap.resize(baseData.numSurfaceNodes,99999.0);
//    surNodes_mSegment.resize(baseData.numSurfaceNodes,-1);
//    surNodes_mSegmentWhichPart.resize(baseData.numSurfaceNodes,-1);
    surNodes_masters.resize(baseData.numSurfaceNodes, 12);
    surNodes_masters.fill(-1);
    segmentsLinkedList.resize(baseData.numSurfaceNodes, 5);
    segmentsLinkedList.fill(-2);
    cellsHeads.resize(numXCells*numYCells, 2);
    cellsHeads.fill(-1);

    
    auto t1 = std::chrono::high_resolution_clock::now();
    
    update_cells(baseData, pars);

    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d1 = t2 - t1;
    std::cout << "elapsed time in verlet cells update:  " << d1.count() << std::endl;
    
    for (int cellYid=0; cellYid<numYCells; cellYid++)
    {
        
        for (int cellXid=0; cellXid<numXCells; cellXid++)
        {
            int cellFlatId = numXCells*cellYid + cellXid + baseData.numSurfaceNodes;
            int slaveNodeId = nodesLinkedList[cellFlatId];

            
            while (slaveNodeId>=0)
            {
                int slaveMesh = baseData.nodeToSegments[baseData.flatSurfaceNodes[slaveNodeId]][2];
                
                for (auto const& delta: neighborBinDelta)
                {
                    int  nCellXid = cellXid + delta.first;
                    int  nCellYid = cellYid + delta.second;
                    
                    if ((nCellXid < 0) || (nCellYid < 0) || (nCellXid == numXCells) || (nCellYid == numYCells))
                    {
                        continue;
                    }
                    
                    int nCellFlatId = numXCells*nCellYid + nCellXid;
                    int masterSegment = cellsHeads(nCellFlatId,0);
                    int nextMSegment;
                    int column = cellsHeads(nCellFlatId,1);
                    
                    while (masterSegment>=0)
                    {
                        int masterMesh = baseData.surfaceSegments[masterSegment][2];
                        
                        if (slaveMesh==masterMesh)
                        {
                            
                            nextMSegment = segmentsLinkedList(masterSegment,column);
                            column = segmentsLinkedList(masterSegment,column+1);
                            masterSegment = nextMSegment;
                            continue;
                        }

                        NTS_interaction(slaveNodeId,masterSegment,masterMesh, baseData, pars);
                        
                        nextMSegment = segmentsLinkedList(masterSegment,column);
                        column = segmentsLinkedList(masterSegment,column+1);
                        masterSegment = nextMSegment;
                    }
                }
                
                slaveNodeId = nodesLinkedList[slaveNodeId];
            }
        }
    }
  
    
    auto t3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d2 = t3 - t2;
    std::cout << "elapsed time in looping over verlet cells:  " << d2.count() << std::endl;
    
    for (int nodeID = 0; nodeID <baseData.numSurfaceNodes; nodeID++){
        
        for(int meshColInc = 0; meshColInc < 3; meshColInc++){
            int meshCol = meshColInc*4;
            if(surNodes_masters(nodeID,meshCol)==-1){
                continue;
            }
            
            int whichPart =  surNodes_masters(nodeID,meshCol+2);
            double gap = surNodes_masters(nodeID,meshCol+3);
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
            int segment = surNodes_masters(nodeID,meshCol+1);
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
                }
                if ( node0 < baseData.numOriginalNodes){
                    forceX(node0) = forceX(node0) + f0x ;
                    forceY(node0) = forceY(node0) + f0y ;
                }
                if ( node1 < baseData.numOriginalNodes){
                    forceX(node1) = forceX(node1) + f1x ;
                    forceY(node1) = forceY(node1) + f1y ;
                }
                contactsEnergy += pars.penaltyStiffness/2 *(gap*gap);
                segmentIinteractions++ ;
                
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
                }
                if ( node0 < baseData.numOriginalNodes){
                    
                    forceX(node0) = forceX(node0) + f00x ;
                    forceY(node0) = forceY(node0) + f00y ;
                }
                contactsEnergy += pars.penaltyStiffness/2 *(gap*gap);
                nodeIinteractions++ ;
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
                }
                if ( node1 < baseData.numOriginalNodes){
                    
                    forceX(node1) = forceX(node1) + f11x ;
                    forceY(node1) = forceY(node1) + f11y ;
                }
                contactsEnergy += pars.penaltyStiffness/2 *(gap*gap);
                nodeIinteractions++ ;
            }
        }
    }
    
    auto t4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d4 = t4 - t3;
    std::cout << "elapsed time in looping over gaps map:  " << d4.count() << std::endl;
}
