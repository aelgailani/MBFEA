#ifndef BaseSysData_hpp
#define BaseSysData_hpp

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


class BaseSysData {
public:
    BaseSysData(const Parameters& pars);
    
    Eigen::Matrix<double,Eigen::Dynamic, 1> refPosX;
    Eigen::Matrix<double,Eigen::Dynamic, 1> refPosY;
    double lyRef;
    double lxRef;
    std::map<int,std::vector<int>> surfaceMeshes;
    std::map<int,std::vector<int>> originalSurfaceMeshes;
    std::vector<int> flatSurfaceNodes;
    std::map<int,std::vector<int>> imageMeshesL;
    std::map<int,std::vector<int>> imageMeshesR;
    std::map<int,std::vector<int>> imageMeshesB;
    std::map<int,std::vector<int>> imageMeshesT;
    std::map<int,std::vector<int>> imageMeshesBL;
    std::map<int,std::vector<int>> imageMeshesBR;
    std::map<int,std::vector<int>> imageMeshesTL;
    std::map<int,std::vector<int>> imageMeshesTR;
    
    std::map<int,std::vector<int>>  triangles;
    std::map<int,std::vector<int>>  augTriangles;
    int numNodes;
    int numOriginalNodes;
    int numSurfaceNodes;
    int numOriginalSurfaceNodes;
    int numMeshes;
    int numOriginalMeshes;
    int numElements;
    int numNodesPerMesh;
    std::vector< std::vector<int> > surfaceSegments;
    std::vector< std::vector<int> > nodeToSegments;
    std::vector<int> nodeToMesh;
    void dump_augmented_surface_meshes(const Parameters& pars);
    
};

#endif /* BaseSysData_hpp */
