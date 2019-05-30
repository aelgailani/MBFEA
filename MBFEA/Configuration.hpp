#ifndef Configuration_hpp
#define Configuration_hpp

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

class Configuration {

public:
    
    Configuration(const BaseSysData& baseData, const Parameters& pars);

    Eigen::VectorXd curPosX;
    Eigen::VectorXd curPosY;
    Eigen::VectorXd curPosXAtLastGridUpdate;
    Eigen::VectorXd curPosYAtLastGridUpdate;
    Eigen::VectorXd augmentedCurPosX;
    Eigen::VectorXd augmentedCurPosY;
    Eigen::VectorXd defGradXX;
    Eigen::VectorXd defGradXY;
    Eigen::VectorXd defGradYX;
    Eigen::VectorXd defGradYY;
    Eigen::VectorXd invDefGradTransXX;
    Eigen::VectorXd invDefGradTransXY;
    Eigen::VectorXd invDefGradTransYX;
    Eigen::VectorXd invDefGradTransYY;
    Eigen::VectorXd forceX;
    Eigen::VectorXd forceY;
    Eigen::VectorXd interForceX;
    Eigen::VectorXd interForceY;
    Eigen::VectorXd PK1stressXX;
    Eigen::VectorXd PK1stressXY;
    Eigen::VectorXd PK1stressYX;
    Eigen::VectorXd PK1stressYY;
    Eigen::VectorXd CstressXX;
    Eigen::VectorXd CstressXY;
    Eigen::VectorXd CstressYX;
    Eigen::VectorXd CstressYY;
    Eigen::VectorXd fSquared;
    Eigen::VectorXd swellingPressurePerEle;
    Eigen::VectorXd elasticEnergyPerEle;
    Eigen::VectorXd mixingEnergyPerEle;
    Eigen::VectorXd internalEnergyPerEle;
    Eigen::VectorXd refArea;
    Eigen::VectorXd areaRatio;
    Eigen::SparseMatrix<double> gradX;
    Eigen::SparseMatrix<double> gradY;
    Eigen::VectorXd wallForceTop;
    Eigen::VectorXd wallForceBottom;
    Eigen::VectorXd wallForceRight;
    Eigen::VectorXd wallForceLeft;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosXL;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosXR ;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosXB;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosXT;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosXBL;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosXBR ;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosXTL;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosXTR;
    
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosYL;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosYR ;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosYB;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosYT ;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosYBL;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosYBR ;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosYTL;
    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosYTR ;
    std::vector<int> masterSlave;
    std::pair<int,int> neighborBinDelta[9] = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,0},{0,1},{1,-1},{1,0},{1,1}};
    double totalEnergy=0, internalEnergy=0, wallsEnergy=0, contactsEnergy=0, shearVirial=0, pressureVirial=0;
    double topPos, botPos, leftPos, rightPos;
    double xMid, yMid;
    double lyCur;
    double lxCur;
    double lyNew;
    double lxNew;
    double A;
    double e0;
    double e1;
    double phi;
    double maxR;
    double LX;
    double LY;
    double Fh;
    double Fv;
    double P1;
    double P2;
    double S1;
    double S2;
    double S3;
    double S4;
    double ex;
    double ey;
    std::map<std::pair<int,int>, std::vector<int>> spatialGridNodes,spatialGridSegments,spatialGridMeshes;
    std::map<std::pair<int,int>, std::vector<double>> gaps;
    std::map<std::pair<int,int>, std::vector<int>> augmentedSegments;
    std::map<std::pair<int,int>, std::vector<int>> augmentedMeshes;
    
    Eigen::VectorXd displacementSinceLastGridUpdate;

    
    void update_post_processing_data(const BaseSysData& baseData, const Parameters& pars);
    void dump_per_node(const BaseSysData& baseData, const Parameters& pars, int& timeStep);
    void dump_per_ele(const BaseSysData& baseData, const Parameters& pars, int& timeStep);
    void compute_forces_walls(const BaseSysData& baseData, const Parameters& pars, const int& timeStep);
    void compute_forces_PBC(const BaseSysData& baseData, const Parameters& pars, const int& timeStep);
    void shear(const BaseSysData& baseData, const Parameters& pars, double strain);
    void compress(const BaseSysData& baseData, const Parameters& pars, double strain);
    void hold(const BaseSysData& baseData, const Parameters& pars);
    void dump_global_data(const Parameters& pars, char mode, char purpose);  //mode: "w" for writing or "a" for appending. purpose: "i" for inspection or "f" for final results
};

#endif /* Configuration_hpp */
