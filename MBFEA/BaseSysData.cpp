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
    sb >> numOriginalNodes;
    refPosX.resize(numOriginalNodes, 1);
    refPosY.resize(numOriginalNodes, 1);
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
    assert(numOriginalNodes==refPosX.size());
    assert(numOriginalNodes==refPosY.size());
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
//    int i =0;
    assert(line=="Number of surfaces");
    std::getline(inFile, line);
    std::istringstream sb1(line);
    sb1 >> numOriginalMeshes;
    std::getline(inFile, line);
    std::istringstream sb2(line);
    sb2 >> numNodesPerMesh;
    numOriginalSurfaceNodes=numNodesPerMesh*numOriginalMeshes;

    while(meshID < numOriginalMeshes)
    {
        while(std::getline(inFile, line))
        {
            if(line==""){break;}
            std::istringstream split(line);
            split >> nodeID;
            originalSurfaceMeshes[meshID].push_back(nodeID);
        }
        meshID++;
    }
    inFile.close();
    
    if (pars.boundaryType == "periodic"){
        for (int meshID=0; meshID < numOriginalMeshes; meshID++){
            for(int& d : originalSurfaceMeshes[meshID] ){
                surfaceMeshes[meshID].push_back(d);
                surfaceMeshes[meshID+numOriginalMeshes].push_back(d+numOriginalNodes);
                surfaceMeshes[meshID+numOriginalMeshes*2].push_back(d+numOriginalNodes*2);
                surfaceMeshes[meshID+numOriginalMeshes*3].push_back(d+numOriginalNodes*3);
                surfaceMeshes[meshID+numOriginalMeshes*4].push_back(d+numOriginalNodes*4);
                surfaceMeshes[meshID+numOriginalMeshes*5].push_back(d+numOriginalNodes*5);
                surfaceMeshes[meshID+numOriginalMeshes*6].push_back(d+numOriginalNodes*6);
                surfaceMeshes[meshID+numOriginalMeshes*7].push_back(d+numOriginalNodes*7);
                surfaceMeshes[meshID+numOriginalMeshes*8].push_back(d+numOriginalNodes*8);
                
            }
        }
        
        numNodes = numOriginalNodes * 9;
        numMeshes = numOriginalMeshes * 9;
        numSurfaceNodes = numOriginalSurfaceNodes * 9;

        for (int meshID=0; meshID < numMeshes; meshID++){
            for(int& d : surfaceMeshes[meshID] ){

                flatSurfaceNodes.push_back(d);
                nodeToMesh.push_back(meshID);
            }
            
            
        }
        
        
        
        
        
        
    }else if (pars.boundaryType == "walls"){
        std::cout << pars.boundaryType<< std::endl;
        numNodes = numOriginalNodes;
        numMeshes = numOriginalMeshes;
        numSurfaceNodes = numOriginalSurfaceNodes;
        surfaceMeshes = originalSurfaceMeshes;
        for (int meshID=0; meshID < numMeshes; meshID++){
            for(int& d : surfaceMeshes[meshID] ){
                
                flatSurfaceNodes.push_back(d);
                nodeToMesh.push_back(meshID);
            }
        }
    }else{
        printf("Please specifiy a valid boundries type ['periodic' or 'walls']!");
        exit(1);
    }
   
    

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
        
        if (pars.boundaryType == "periodic"){
            augTriangles[eleID].push_back(n1);
            augTriangles[eleID].push_back(n2);
            augTriangles[eleID].push_back(n3);
            
            augTriangles[eleID+numElements].push_back(n1+numOriginalNodes);
            augTriangles[eleID+numElements].push_back(n2+numOriginalNodes);
            augTriangles[eleID+numElements].push_back(n3+numOriginalNodes);
            augTriangles[eleID+numElements*2].push_back(n1+numOriginalNodes*2);
            augTriangles[eleID+numElements*2].push_back(n2+numOriginalNodes*2);
            augTriangles[eleID+numElements*2].push_back(n3+numOriginalNodes*2);
            augTriangles[eleID+numElements*3].push_back(n1+numOriginalNodes*3);
            augTriangles[eleID+numElements*3].push_back(n2+numOriginalNodes*3);
            augTriangles[eleID+numElements*3].push_back(n3+numOriginalNodes*3);
            augTriangles[eleID+numElements*4].push_back(n1+numOriginalNodes*4);
            augTriangles[eleID+numElements*4].push_back(n2+numOriginalNodes*4);
            augTriangles[eleID+numElements*4].push_back(n3+numOriginalNodes*4);
            augTriangles[eleID+numElements*5].push_back(n1+numOriginalNodes*5);
            augTriangles[eleID+numElements*5].push_back(n2+numOriginalNodes*5);
            augTriangles[eleID+numElements*5].push_back(n3+numOriginalNodes*5);
            augTriangles[eleID+numElements*6].push_back(n1+numOriginalNodes*6);
            augTriangles[eleID+numElements*6].push_back(n2+numOriginalNodes*6);
            augTriangles[eleID+numElements*6].push_back(n3+numOriginalNodes*6);
            augTriangles[eleID+numElements*7].push_back(n1+numOriginalNodes*7);
            augTriangles[eleID+numElements*7].push_back(n2+numOriginalNodes*7);
            augTriangles[eleID+numElements*7].push_back(n3+numOriginalNodes*7);
            augTriangles[eleID+numElements*8].push_back(n1+numOriginalNodes*8);
            augTriangles[eleID+numElements*8].push_back(n2+numOriginalNodes*8);
            augTriangles[eleID+numElements*8].push_back(n3+numOriginalNodes*8);

        }
        
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
    
    if (pars.boundaryType=="periodic"){
        std::ofstream augTrifile;
           augTrifile.open (pars.outputFolderName+"/augElements.txt");
           for (int i=0; i< augTriangles.size();i++ ) {
               augTrifile << augTriangles[i][0] << std::setw(10)<< augTriangles[i][1] << std::setw(10)<< augTriangles[i][2] <<std::endl;
           }
           augTrifile.close();
    }

//******* fill in segments data ********************************************************************************************************************
   
    /* Create segments and there relations to the nodes ** specificly for the input forms in our project **
     surfaceSegments[segmentID] = {firstNodeID, secondNodeID, meshID, prevSegid, nextSegID}
     nodeToSegments[nodeID] = {firstNodeID, secondNodeID, meshID, prevSegid, nextSegID}
     */
    if (pars.contactMethod=="nts" || pars.contactMethod=="gntn" || pars.contactMethod=="gntn2"  ){
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
            surfaceSegments[segmentID-numNodesPerMesh+1][3] = segmentID;
            nodeToSegments[surfaceMeshes[meshID][surfaceNodeID]][1] = segmentID;
            nodeToSegments[surfaceMeshes[meshID][0]][0] = segmentID;
            nodeToSegments[surfaceMeshes[meshID][surfaceNodeID]][2] = meshID;
            segmentID++;
        }
    }





}

void BaseSysData::dump_augmented_surface_meshes(const Parameters& pars){
    std::ofstream myfile;
    myfile.open (pars.outputFolderName+"/augmentedSurfaceMeshes.txt");
    myfile << "Number of surfaces" << "\n" << numMeshes << "\n" << numNodesPerMesh << std::endl;
    for (int meshID=0; meshID < numMeshes;meshID++ ) {
        for(int& d : surfaceMeshes[meshID] ){
            myfile << d << std::endl;
        }
        myfile << "\n";
        
    }
    myfile.close();
}
