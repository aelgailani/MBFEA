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
#include "Parameters.hpp"
#include "BaseSysData.hpp"

BaseSysData::BaseSysData(const Parameters& pars){
    
    lyRef = pars.initTopPos - pars.initBotPos;
    lxRef = pars.initRightPos - pars.initLeftPos;
//**************  Read nodes data  ********************************************************************************************
    std::ifstream inFile;
    inFile.open(pars.initialNodesFileName);
        ///Check for errors
    if (inFile.fail())
    {
        std::cerr << "Error Openning Initial Nodes File" << std::endl;
        exit(1);
    }
        ///Read the number of nodes found in the first line of the input file
    std::string line;
    std::getline(inFile, line);
    std::istringstream sb(line);
    sb >> numNodes;
    refPosX.resize(numNodes, 1);
    refPosY.resize(numNodes, 1);
    double ID,isBoundary,x,y;
    int id = 0;
    while(std::getline(inFile, line))
    {
        std::istringstream split(line);
        split >> ID >> x >> y >> isBoundary;
        refPosX(id) = x;
        refPosY(id) = y;
        id++;
    }
    assert(numNodes==refPosX.size());
    assert(numNodes==refPosY.size());
    inFile.close();

    
//**************  Read Surface nodes data  ********************************************************************************************

    inFile.open(pars.surfaceNodesFileName);
    //Check for errors
    if (inFile.fail())
    {
        std::cerr << "Error Openning Surface Nodes File" << std::endl;
        exit(1);
    }
    //Read the number of meshes found in the first line of the input file
    std::getline(inFile, line);
    int meshID = 0;
    int nodeID;
    int i =0;
    assert(line=="Number of surfaces");
    std::getline(inFile, line);
    std::istringstream sb1(line);
    sb1 >> numMeshes;
    std::getline(inFile, line);
    std::istringstream sb2(line);
    sb2 >> numNodesPerMesh;
    
    while(meshID < numMeshes)
    {
        i =0;
        while(std::getline(inFile, line))
        {
            if(line==""){break;}
            std::istringstream split(line);
            split >> nodeID;
            surfaceMeshes[meshID].push_back(nodeID);
            i++;
            (numSurfaceNodes)++;
        }
        meshID++;
    }
    inFile.close();
    
    
//**************  Read triangles data  ********************************************************************************************
    
    inFile.open(pars.trianglesFileName); //Make sure later that the file is open in "read" mode only
    //Check for errors
    if (inFile.fail())
    {
        std::cerr << "Error Openning Elements File" << std::endl;
        exit(1);
    }
    //Read the number of elements found in the first line of the input file
    std::getline(inFile, line);
    std::istringstream sb3(line);
    sb3 >> numElements;
    int eleID = 0;
    int n1,n2,n3;
    while(eleID < numElements)
    {
        std::getline(inFile, line);
        if(line==""){break;}
        std::istringstream split(line);
        split >> nodeID >> n1 >> n2 >> n3;
        triangles[eleID].push_back(n1);
        triangles[eleID].push_back(n2);
        triangles[eleID].push_back(n3);
        eleID++;
    }
    inFile.close();

//**************  Dump triangels to the output folder ********************************************************************************************
    std::ofstream trifile;
    trifile.open (pars.outputFolderName+"/elements.txt");
    for (int i=0; i< triangles.size();i++ ) {
        trifile << triangles[i][0] << std::setw(10)<< triangles[i][1] << std::setw(10)<< triangles[i][2] <<std::endl;
    }
    trifile.close();

    
//******* fill in segments data ********************************************************************************************************************
   
    /* Create segments and there relations to the nodes ** specificly for the input forms in our project **
     surfaceSegments[segmentID] = {firstNodeID, secondNodeID, meshID, prevSegid, nextSegID}
     nodeToSegments[nodeID] = {firstNodeID, secondNodeID, meshID, prevSegid, nextSegID}
     */
    
    surfaceSegments.resize(numSurfaceNodes,std::vector<int>(5));
    nodeToSegments.resize(numNodes,std::vector<int>(3,999999999));
    int segmentID = 0;
    int surfaceNodeID = 0;
    for (int meshID = 0; meshID < numMeshes; meshID++)
    {
        surfaceNodeID = 0;
        for (surfaceNodeID = 0; surfaceNodeID < (numNodesPerMesh-1) ; surfaceNodeID++)
        {
            surfaceSegments[segmentID][0] = surfaceMeshes[meshID][surfaceNodeID];
            surfaceSegments[segmentID][1] = surfaceMeshes[meshID][surfaceNodeID + 1];
            surfaceSegments[segmentID][2] = meshID;
            surfaceSegments[segmentID][3] = segmentID-1;
            surfaceSegments[segmentID][4] = segmentID+1;
            nodeToSegments[surfaceMeshes[meshID][surfaceNodeID+1]][0] = segmentID;
            nodeToSegments[surfaceMeshes[meshID][surfaceNodeID]][1] = segmentID;
            nodeToSegments[surfaceMeshes[meshID][surfaceNodeID]][2] = meshID;
            segmentID++;
        }
        surfaceSegments[segmentID][0] = surfaceMeshes[meshID][surfaceNodeID];
        surfaceSegments[segmentID][1] = surfaceMeshes[meshID][0];
        surfaceSegments[segmentID][2] = meshID;
        surfaceSegments[segmentID][3] = segmentID-1;
        surfaceSegments[segmentID][4] = segmentID-numNodesPerMesh+1;
        surfaceSegments[0][3] = segmentID; // check later. looks weird.
        nodeToSegments[surfaceMeshes[meshID][surfaceNodeID]][1] = segmentID;
        nodeToSegments[surfaceMeshes[meshID][0]][0] = segmentID;
        nodeToSegments[surfaceMeshes[meshID][surfaceNodeID]][2] = meshID;
        segmentID++;
    }





}

