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
    
    numXBins = floor(lxNew*(1+2*pars.imagesMargin)/pars.verletCellCutoff);
    numYBins = floor(lyNew*(1+2*pars.imagesMargin)/pars.verletCellCutoff);
    verletCellSizeX  = lxNew*(1+2*pars.imagesMargin)/numXBins;
    verletCellSizeY = lyNew*(1+2*pars.imagesMargin)/numYBins;
    cellListNodes.resize(numXBins*numYBins+baseData.numSurfaceNodes,-1);

    cellListSegments.resize(baseData.numSurfaceNodes+1, numXBins*numYBins);
    cellListSegments.fill(-2);
    
    for (int i=0; i < numXBins*numYBins; i++)
    {
        cellListSegments(baseData.numSurfaceNodes,i)=-1;
    }
    
    auto t1 = std::chrono::high_resolution_clock::now();
    
    update_cells(baseData, pars);
    
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d1 = t2 - t1;
    std::cout << "elapsed time in verlet cells update:  " << d1.count() << std::endl;
    
    for (int m1y=0; m1y<numYBins; m1y++)
    {
        
        for (int m1x=0; m1x<numXBins; m1x++)
        {
            int m1 = numXBins*m1y + m1x + baseData.numSurfaceNodes;
            
            int slaveNodeId = cellListNodes[m1];
            int slaveMesh = baseData.nodeToSegments[baseData.flatSurfaceNodes[slaveNodeId]][2];
            while (slaveNodeId>=0)
            {
                for (auto const& delta: neighborBinDelta)
                {
                    int  nBinXid = m1x + delta.first;
                    int  nBinYid = m1y + delta.second;
                    
                    if ((nBinXid < 0) || (nBinYid < 0) || (nBinXid == numXBins) || (nBinYid == numYBins))
                    {
                        continue;
                        
                    }
                    int   m2 = numXBins*nBinYid + nBinXid;
                    int masterSegment = cellListSegments(baseData.numSurfaceNodes,m2);
                    while (masterSegment>=0)
                    {
                        int masterMesh = baseData.surfaceSegments[masterSegment][2];
                        if (slaveMesh==masterMesh)
                        {
                            masterSegment = cellListSegments(masterSegment,m2);
                            continue;
                        }

                        NTS_interaction(baseData.flatSurfaceNodes[slaveNodeId],masterSegment, baseData, pars);
                        
                        masterSegment = cellListSegments(masterSegment,m2);
                    }
                }
                
                slaveNodeId = cellListNodes[slaveNodeId];
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
