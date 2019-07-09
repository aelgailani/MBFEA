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
    gaps.clear();
    
    numXCells = floor(lxNew*(1+2*pars.imagesMargin)/pars.verletCellCutoff);
    numYCells = floor(lyNew*(1+2*pars.imagesMargin)/pars.verletCellCutoff);
    verletCellSizeX  = lxNew*(1+2*pars.imagesMargin)/numXCells;
    verletCellSizeY = lyNew*(1+2*pars.imagesMargin)/numYCells;
    
    nodesLinkedList.resize(numXCells*numYCells+baseData.numSurfaceNodes,-1);
    segmentsLinkedList.resize(baseData.numSurfaceNodes, std::vector<int>(5, -2));
    cellsHeads.resize(numXCells*numYCells, std::vector<int>(2, -1));
    
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
                    int masterSegment = cellsHeads[nCellFlatId][0];
                    int nextMSegment;
                    int column = cellsHeads[nCellFlatId][1];
                    
                    while (masterSegment>=0)
                    { 
                        int masterMesh = baseData.surfaceSegments[masterSegment][2];
                        
                        if (slaveMesh==masterMesh)
                        {
                            
                            nextMSegment = segmentsLinkedList[masterSegment][column];
                            column = segmentsLinkedList[masterSegment][column+1];
                            masterSegment = nextMSegment;
                            continue;
                        }

                        NTS_interaction(baseData.flatSurfaceNodes[slaveNodeId],masterSegment, baseData, pars);
                        
                        nextMSegment = segmentsLinkedList[masterSegment][column];
                        column = segmentsLinkedList[masterSegment][column+1];
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
    
    for (auto const& nodeRow: gaps)
    {
        
        std::vector<double> shortestPath = nodeRow.second;
        
        if (shortestPath[1]>= 0)
        {
            continue;
        }
        if (shortestPath[0] > maxInterference)
        {
            maxInterference = shortestPath[0];
        }
        if (shortestPath.size()==20)
        {
            
            segmentIinteractions++ ;
            if (nodeRow.first.first < baseData.numOriginalNodes)
            {
                forceX(nodeRow.first.first) = forceX(nodeRow.first.first) + shortestPath[4] ;
                forceY(nodeRow.first.first) = forceY(nodeRow.first.first) + shortestPath[5] ;
            }
            if ( shortestPath[6] < baseData.numOriginalNodes){
                forceX(shortestPath[6]) = forceX(shortestPath[6])  + shortestPath[8] ;
                forceY(shortestPath[6]) = forceY(shortestPath[6]) + shortestPath[9];
            }
            if ( shortestPath[10] < baseData.numOriginalNodes){
                forceX(shortestPath[10]) = forceX(shortestPath[10]) + shortestPath[12];
                forceY(shortestPath[10]) = forceY(shortestPath[10]) + shortestPath[13];
            }
            
            
            contactsEnergy += pars.penaltyStiffness/2 *(shortestPath[0]*shortestPath[0]);
            
            
        }else{
            
            nodeIinteractions++;
            if (nodeRow.first.first < baseData.numOriginalNodes){
                forceX(nodeRow.first.first) = forceX(nodeRow.first.first) + shortestPath[4] ;
                forceY(nodeRow.first.first) = forceY(nodeRow.first.first) + shortestPath[5] ;
            }
            if ( shortestPath[6] < baseData.numOriginalNodes){
                forceX(shortestPath[6]) = forceX(shortestPath[6])  + shortestPath[8] ;
                forceY(shortestPath[6]) = forceY(shortestPath[6]) + shortestPath[9];
            }
            
            contactsEnergy += pars.penaltyStiffness/2 *(shortestPath[0]*shortestPath[0]);
            
        }
    }
    
    auto t4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d4 = t4 - t3;
    std::cout << "elapsed time in looping over gaps map:  " << d4.count() << std::endl;
}
