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
#include <valarray>
#include "Parameters.hpp"
#include "BaseSysData.hpp"

class Configuration {

public:
    
    Configuration(const BaseSysData& baseData, const Parameters& pars);

    Eigen::VectorXd curPosX;
    Eigen::VectorXd curPosY;
    Eigen::VectorXd curPosXAtLastStep;
    Eigen::VectorXd curPosYAtLastStep;
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
    Eigen::VectorXd velocityX;
    Eigen::VectorXd velocityY;
    Eigen::VectorXd prevVelocityX;
    Eigen::VectorXd prevVelocityY;
    Eigen::VectorXd interForceX;
    Eigen::VectorXd interForceY;
    Eigen::VectorXd surfaceForceX;
    Eigen::VectorXd surfaceForceY;
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
    
    Eigen::VectorXd curPosXL;
    Eigen::VectorXd curPosXR ;
    Eigen::VectorXd curPosXB;
    Eigen::VectorXd curPosXT;
    Eigen::VectorXd curPosXBL;
    Eigen::VectorXd curPosXBR ;
    Eigen::VectorXd curPosXTL;
    Eigen::VectorXd curPosXTR;
    
    Eigen::VectorXd curPosYL;
    Eigen::VectorXd curPosYR ;
    Eigen::VectorXd curPosYB;
    Eigen::VectorXd curPosYT ;
    Eigen::VectorXd curPosYBL;
    Eigen::VectorXd curPosYBR ;
    Eigen::VectorXd curPosYTL;
    Eigen::VectorXd curPosYTR;
    
    Eigen::VectorXd consistencyFactorX;
    Eigen::VectorXd consistencyFactorY;
    Eigen::VectorXd consistencyErrorFactorX;
    Eigen::VectorXd consistencyErrorFactorY;
    
    int segmentIinteractions, nodeIinteractions;
    double totalEnergy=0, internalEnergy=0, wallsEnergy=0, contactsEnergy=0, shearVirial=0, pressureVirial=0;
    double prevTotEnergy=0, deltaTotEnergy=0, L2NormResidual=0 ;
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
    double avgR;
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
    int numXCells, numYCells;
    double verletCellSizeX;
    double verletCellSizeY;
    double maxWallinterference;
    double maxInterference;
//    abs(gap),gapSign,double(node),f,fx,fy,double(node0),f0,f0x,f0y,double(node1),f1,f1x,f1y,nx,ny,s,double(segment),(xi-x0)*nx,(yi-y0)*ny
    std::map<std::pair<int,int>, std::vector<double>> gaps;
    std::valarray<double> surNodes_gap;
    std::valarray<int> surNodes_mSegment;
    std::valarray<int> surNodes_mSegmentWhichPart;
    std::valarray<int> nodesLinkedList;

    std::valarray<int> surNodes_mMesh1;
    std::valarray<int> surNodes_mMesh2;
    std::valarray<int> surNodes_mMesh3;
    std::valarray<int> surNodes_mSegment1;
    std::valarray<int> surNodes_mSegment2;
    std::valarray<int> surNodes_mSegment3;
    std::valarray<int> surNodes_mPart1;
    std::valarray<int> surNodes_mPart2;
    std::valarray<int> surNodes_mPart3;
    std::valarray<double> surNodes_gap1;
    std::valarray<double> surNodes_gap2;
    std::valarray<double> surNodes_gap3;
    
    Eigen::MatrixXd segmentsLinkedList;
    Eigen::MatrixXd cellsHeads;
    
    std::vector<int> masterSlave;
    std::pair<int,int> neighborBinDelta[9] = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,0},{0,1},{1,-1},{1,0},{1,1}};
    
    Eigen::VectorXd displacementSinceLastStep;
    void check_force_energy_consistency(const BaseSysData& baseData, const Parameters& pars);
    void update_cells(const BaseSysData& baseData, const Parameters& pars);
    void update_post_processing_data(const BaseSysData& baseData, const Parameters& pars);
    void dump_per_node(const BaseSysData& baseData, const Parameters& pars, int& timeStep);
    void dump_per_node_periodic_images_on(const BaseSysData& baseData, const Parameters& pars, int& timeStep);
    void dump_per_ele(const BaseSysData& baseData, const Parameters& pars, int& timeStep);
    void compute_forces_walls(const BaseSysData& baseData, const Parameters& pars, const int& timeStep);
    void compute_forces_PBC(const BaseSysData& baseData, const Parameters& pars, const int& timeStep, bool surfaceInteractions, bool updatePBC);
    void compute_surface_forces(const BaseSysData& baseData, const Parameters& pars, const int& timeStep);
    void NTS_interaction(const int& node, const int& segment,const int& masterMesh, const BaseSysData& baseData, const Parameters& pars);
    void shear(const BaseSysData& baseData, const Parameters& pars, double strain);
    void special_localized_deformation(const BaseSysData& baseData, const Parameters& pars,const double& gammaX, const double& gammaY, const std::vector<int>& targetNodes);
    void compress(const BaseSysData& baseData, const Parameters& pars, double strain);
    void hold(const BaseSysData& baseData, const Parameters& pars);
    void dump_global_data(const Parameters& pars, char mode, char purpose);  //mode: "w" for writing or "a" for appending. purpose: "i" for inspection or "f" for final results
};

#endif /* Configuration_hpp */
