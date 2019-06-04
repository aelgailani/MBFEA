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

    
    numXBins = floor(lxNew*(1+2*pars.imagesMargin)/pars.verletCellCutoff);
    numYBins = floor(lyNew*(1+2*pars.imagesMargin)/pars.verletCellCutoff);
    verletCellSizeX  = lxNew*(1+2*pars.imagesMargin)/numXBins;
    verletCellSizeY = lyNew*(1+2*pars.imagesMargin)/numYBins;
    cellList.resize(numXBins*numYBins+baseData.numSurfaceNodes,-1);
    
    auto t1 = std::chrono::high_resolution_clock::now();
    
    update_cells(baseData, pars);
    
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d1 = t2 - t1;
    std::cout << "elapsed time in verlet cells update:  " << d1.count() << std::endl;
    std::chrono::duration<double> d222 = t2-t2;
    
    for (int m1y=0; m1y<numYBins; m1y++)
    {
        
        for (int m1x=0; m1x<numXBins; m1x++)
        {
            int m1 = numXBins*m1y + m1x + baseData.numSurfaceNodes;
            
            closestMaster.resize(11,9999);
            
            for (auto const& delta: neighborBinDelta)
            {
                int  nBinXid = m1x + delta.first;
                int  nBinYid = m1y + delta.second;
                
                if ((nBinXid < 0) || (nBinYid < 0) || (nBinXid == numXBins) || (nBinYid == numYBins))
                {
                    continue;
                    
                }
                int   m2 = numXBins*nBinYid + nBinXid + baseData.numSurfaceNodes;
                int slaveNodeId = cellList[m1];  //this is the slave surface node id in flatSurfaceNodes vector, not the global id
                int slaveNode = baseData.flatSurfaceNodes[slaveNodeId]; //this is the slave surface node id globally
                while (slaveNodeId>=0)
                {
                    int slaveMesh = baseData.nodeToSegments[baseData.flatSurfaceNodes[slaveNodeId]][2];
                    int masterNodeId = cellList[m2]; //this is the master surface node id in flatSurfaceNodes vector, not the global id
                    while (masterNodeId>=0)
                    {
                        int masterMesh = baseData.nodeToSegments[baseData.flatSurfaceNodes[masterNodeId]][2];
                        if (slaveMesh==masterMesh)
                        {
                            masterNodeId = cellList[masterNodeId];
                            continue;
                        }
                        
                        auto t222 = std::chrono::high_resolution_clock::now();
                        
                        int segment0 = baseData.nodeToSegments[baseData.flatSurfaceNodes[masterNodeId]][0];
                        int segment1 = baseData.nodeToSegments[baseData.flatSurfaceNodes[masterNodeId]][1];
                        NTS_interaction(slaveNode,segment0, baseData, pars);
                        NTS_interaction(slaveNode,segment1, baseData, pars);
                        
                        auto t333 = std::chrono::high_resolution_clock::now();
                        std::chrono::duration<double> d2222 = t333 - t222;
                        d222+=d2222;
                        
                        
                        masterNodeId = cellList[masterNodeId];
                    }
                    
                    if (closestMaster[1] >= 0)
                    {
                        closestMaster.resize(11,9999);
                        slaveNodeId = cellList[slaveNodeId];
                        continue;
                    }
                    
                    if (slaveNode< baseData.numOriginalNodes)
                    {
                        forceX(slaveNode) = forceX(slaveNode) + closestMaster[2] ;
                        forceY(slaveNode) = forceY(slaveNode) + closestMaster[3] ;
                    }
                    if ( closestMaster[4] < baseData.numOriginalNodes){
                        forceX(closestMaster[4]) = forceX(closestMaster[4])  + closestMaster[5] ;
                        forceY(closestMaster[4]) = forceY(closestMaster[4]) + closestMaster[6];
                    }
                    if ( closestMaster[7] < baseData.numOriginalNodes){
                        forceX(closestMaster[7]) = forceX(closestMaster[7]) + closestMaster[8];
                        forceY(closestMaster[7]) = forceY(closestMaster[7]) + closestMaster[9];
                    }
                    
                    
                    contactsEnergy += pars.penaltyStiffness/2 *(closestMaster[0]*closestMaster[0]);
                    
                    if (closestMaster[10]==1)
                    {
                        segmentIinteractions++;
                    }else{
                        nodeIinteractions++;
                    }
                    closestMaster.resize(11,9999);
                    slaveNodeId = cellList[slaveNodeId];
                }
            }
        }
    }
  
    auto t3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d2 = t3 - t2;
    std::cout << "elapsed time in looping over verlet cells:  " << d2.count() << std::endl;
    std::cout << "elapsed time in doing NTS arithmatics:  " << d222.count() << std::endl;
    
    
    auto t4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d4 = t4 - t3;
    std::cout << "elapsed time in looping over gaps map:  " << d4.count() << std::endl;
}
