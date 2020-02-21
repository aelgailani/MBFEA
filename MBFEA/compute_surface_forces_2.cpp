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

void Configuration::compute_surface_forces_2(const BaseSysData& baseData, const Parameters& pars, const int& timeStep)
{
    maxInterference = 0;
    segmentIinteractions = 0;
    nodeIinteractions = 0;

    
    numXCells = floor(lxNew*(1+2*pars.imagesMargin)/pars.verletCellCutoff);
    numYCells = floor(lyNew*(1+2*pars.imagesMargin)/pars.verletCellCutoff);
    verletCellSizeX  = lxNew*(1+2*pars.imagesMargin)/numXCells;
    verletCellSizeY = lyNew*(1+2*pars.imagesMargin)/numYCells;
    
    nodesLinkedList.resize(numXCells*numYCells+baseData.numSurfaceNodes,-1);
    
    surNodes_mMesh1.resize(baseData.numSurfaceNodes,-1);
    surNodes_mMesh2.resize(baseData.numSurfaceNodes,-1);
    surNodes_mMesh3.resize(baseData.numSurfaceNodes,-1);
    surNodes_mSegment1.resize(baseData.numSurfaceNodes,-1);
    surNodes_mSegment2.resize(baseData.numSurfaceNodes,-1);
    surNodes_mSegment3.resize(baseData.numSurfaceNodes,-1);
    surNodes_mPart1.resize(baseData.numSurfaceNodes,-1);
    surNodes_mPart2.resize(baseData.numSurfaceNodes,-1);
    surNodes_mPart3.resize(baseData.numSurfaceNodes,-1);
    surNodes_gap1.resize(baseData.numSurfaceNodes,-1);
    surNodes_gap2.resize(baseData.numSurfaceNodes,-1);
    surNodes_gap3.resize(baseData.numSurfaceNodes,-1);
    
    segmentsLinkedList_2.resize(baseData.numSurfaceNodes, 5);
    segmentsLinkedList_2.fill(-2);
    cellsHeads.resize(numXCells*numYCells, 2);
    cellsHeads.fill(-1);

    
    auto t1 = std::chrono::high_resolution_clock::now();
    
    update_cells_2(baseData, pars);

    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d1 = t2 - t1;
    std::cout << "elapsed time in verlet cells update:  " << d1.count() << std::endl;
    
    for (int cellYid=0; cellYid<numYCells; cellYid++)
    {
        
        for (int cellXid=0; cellXid<numXCells; cellXid++)
        {
            int cellFlatId = numXCells*cellYid + cellXid + baseData.numSurfaceNodes;
            int slaveNodeId = nodesLinkedList[cellFlatId];  // Note that slaveNodeId is its local ordinal number not the golbal node name assigned in the mesh

            
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
                            
                            nextMSegment = segmentsLinkedList_2(masterSegment,column);
                            column = segmentsLinkedList_2(masterSegment,column+1);
                            masterSegment = nextMSegment;
                            continue;
                        }

                        NTS_interaction(slaveNodeId,masterSegment,masterMesh, baseData, pars);
                        
                        nextMSegment = segmentsLinkedList_2(masterSegment,column);
                        column = segmentsLinkedList_2(masterSegment,column+1);
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
        
            if(surNodes_mMesh1[nodeID]==-1){
                continue;
            }
            
            int whichPart =  surNodes_mPart1[nodeID];
            double gap = surNodes_gap1[nodeID];
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
            int segment = surNodes_mSegment1[nodeID];
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
//                    std::cout << "sNode:  " << node << std::endl;
//                    std::cout << "mNode0:  " << node0 << std::endl;
//                    std::cout << "mNode1:  " << node1 << std::endl;
//                    std::cout << "gap:  " << gap << std::endl;
//                    std::cout << "f1ix:  " << f1ix << std::endl;
//                    std::cout << "f1iy:  " << f1iy << std::endl;
//                    std::cout << "f11x:  " << f11x << std::endl;
//                    std::cout << "f11y:  " << f11y << std::endl;
//                    std::cout << "contactsEnergy +=  " << pars.penaltyStiffness/2 *(gap*gap) << std::endl;
//                    std::cout << "contactsEnergy =  " << contactsEnergy << std::endl;
                }
            }
        
    }
    
    for (int nodeID = 0; nodeID <baseData.numSurfaceNodes; nodeID++){
        
        if(surNodes_mMesh2[nodeID]==-1){
            continue;
        }
        
        int whichPart =  surNodes_mPart2[nodeID];
        double gap = surNodes_gap2[nodeID];
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
        int segment = surNodes_mSegment2[nodeID];
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
        }
        
    }
    
    
    for (int nodeID = 0; nodeID <baseData.numSurfaceNodes; nodeID++){
        
        if(surNodes_mMesh3[nodeID]==-1){
            continue;
        }
        
        int whichPart =  surNodes_mPart3[nodeID];
        double gap = surNodes_gap3[nodeID];
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
        int segment = surNodes_mSegment3[nodeID];
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
        }
        
    }
    
    auto t4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> d4 = t4 - t3;
    std::cout << "elapsed time in looping over gaps map:  " << d4.count() << std::endl;
}
