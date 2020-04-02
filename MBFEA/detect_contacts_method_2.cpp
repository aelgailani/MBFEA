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

void Configuration::detect_contacts_method_2(const BaseSysData& baseData, const Parameters& pars)
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

                          nts_interaction(slaveNodeId,masterSegment,masterMesh, baseData, pars);
                          
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
