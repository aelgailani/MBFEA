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
#include "utility_functions.hpp"

void Configuration::detect_nts_contacts_two_points_method(const BaseSysData& baseData, const Parameters& pars)
{
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

                           if (not pars.reversibleMasterSlaveRole)
                             {
                                 // slaveMaster is 1 if the slaveMesh is enslaved to the masterMesh, -2 if it is actually the master, and -1 if no claim is etablished yet
                                 if (slaveMaster(slaveMesh,masterMesh)!=1){
                                     
                                     if (slaveMaster(slaveMesh,masterMesh)==-2){
                                             nextMSegment = segmentsLinkedList_2(masterSegment,column);
                                             column = segmentsLinkedList_2(masterSegment,column+1);
                                             masterSegment = nextMSegment;
                                             continue;
                                         }else if (slaveMaster(slaveMesh,masterMesh)==-1){
                                             slaveMaster(slaveMesh,masterMesh)=1;
                                             slaveMaster(masterMesh,slaveMesh)=-2;
                                         }
                                 }

                             }
                          
                          if (pars.smoothCorners==true){
                              nts_find_closest_approach_with_smoothing(slaveNodeId,masterSegment,masterMesh, baseData, pars);
                          }else{
                              nts_find_closest_approach(slaveNodeId,masterSegment,masterMesh, baseData, pars);
                          }
                          
                          nextMSegment = segmentsLinkedList_2(masterSegment,column);
                          column = segmentsLinkedList_2(masterSegment,column+1);
                          masterSegment = nextMSegment;
                      }
                  }
                  
                  slaveNodeId = nodesLinkedList[slaveNodeId];
              }
          }
      }
    
}
