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
    std::map<int,std::vector<int>>  triangles;
    int numNodes;
    int numSurfaceNodes;
    int numMeshes;
    int numElements;
    int numNodesPerMesh;
    std::vector< std::vector<int> > surfaceSegments;
    std::vector< std::vector<int> > nodeToSegments;
    
};

#endif /* BaseSysData_hpp */
