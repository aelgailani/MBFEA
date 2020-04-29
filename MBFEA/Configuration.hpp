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
#include <set>

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
    Eigen::VectorXd DVxDx;
    Eigen::VectorXd DVxDy;
    Eigen::VectorXd DVyDx;
    Eigen::VectorXd DVyDy;
    
    
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
    double KWoodXX, KWoodXY, KWoodYX, KWoodYY;
    double xMid, yMid;
    double lyCur;
    double lxCur;
    double lyNew;
    double lxNew;
    double A;
    double e0=0;
    double prev_e0=0;
    double e1=0;
    double prev_e1=0;
    double phi;
    double maxR;
    double avgR;
    double LX;
    double LY;
    double Fh;
    double Fv;
    double P1;
    double P2;
    double prev_P2=0;
    double S1;
    double S2;
    double prev_S2=0;
    double S3;
    double S4;
    double ex;
    double ey;
    int numXCells, numYCells;
    double verletCellSizeX;
    double verletCellSizeY;
    double maxWallinterference;
    double maxInterference;
    double DPOverDe0;
    double DSOverDe1;
    double A_material;
//    abs(gap),gapSign,double(node),f,fx,fy,double(node0),f0,f0x,f0y,double(node1),f1,f1x,f1y,nx,ny,s,double(segment),(xi-x0)*nx,(yi-y0)*ny
//    std::map<std::pair<int,int>, std::vector<double>> gaps;
//    std::valarray<double> surNodes_gap;
//    std::valarray<int> surNodes_mSegment;
//    std::valarray<int> surNodes_mSegmentWhichPart;
    std::valarray<int> nodesLinkedList;
    std::valarray<int> segmentsLinkedList_1;
    Eigen::MatrixXd segmentsLinkedList_2;
    Eigen::MatrixXd cellsHeads;
    std::map<std::pair<int,int>, std::vector<int>> facets;
    Eigen::MatrixXd slaveMaster;//this map gives you infromation about (master, slave) nodes with duplication allowed. notice that key (1,3) is different from (3,1), the former contains all nodes in mesh 1 while the latter the nodes of mesh 3.
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
    
    
    
    std::pair<int,int> neighborBinDelta[9] = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,0},{0,1},{1,-1},{1,0},{1,1}};
    
    Eigen::VectorXd displacementSinceLastStep;
    std::string lastStepFirst ; //for managing output global data.txt
    void check_force_energy_consistency(const BaseSysData& baseData, const Parameters& pars);
    void update_cells_1(const BaseSysData& baseData, const Parameters& pars); // used with segmentCellList method 1
    void update_cells_2(const BaseSysData& baseData, const Parameters& pars);// used with segmentCellList method 2
    void update_post_processing_data(const BaseSysData& baseData, const Parameters& pars);
    void dump_per_node(const BaseSysData& baseData, const Parameters& pars, long& timeStep);
    void dump_per_node_periodic_images_on(const BaseSysData& baseData, const Parameters& pars, long& timeStep);
    void dump_per_ele(const BaseSysData& baseData, const Parameters& pars, long& timeStep);
    void dump_facets(const BaseSysData& baseData, const Parameters& pars, long& timeStep);
    
    void compute_forces_harmonic_walls(const BaseSysData& baseData, const Parameters& pars, const long& timeStep, bool surfaceInteractions, bool updatePBC, bool Hessian);
    void compute_forces_pbc(const BaseSysData& baseData, const Parameters& pars, const long& timeStep, bool surfaceInteractions, bool updatePBC, bool Hessian);
    void compute_surface_forces(const BaseSysData& baseData, const Parameters& pars, bool Hessian, const long& timeStep);
    void detect_nts_contacts_single_point_method(const BaseSysData& baseData, const Parameters& pars);
    void detect_nts_contacts_two_points_method(const BaseSysData& baseData, const Parameters& pars);
    void apply_nts_harmonic_penalty(const BaseSysData& baseData, const Parameters& pars, const std::valarray<int>& surNodes_mMesh, const std::valarray<int>& surNodes_mSegment, const std::valarray<int>& surNodes_mPart, const std::valarray<double>& surNodes_gap, bool Hessian, const long& timeStep);
    void apply_nts_powerlaw_penalty(const BaseSysData& baseData, const Parameters& pars, const std::valarray<int>& surNodes_mMesh, const std::valarray<int>& surNodes_mSegment, const std::valarray<int>& surNodes_mPart, const std::valarray<double>& surNodes_gap, bool Hessian, const long& timeStep);

    void apply_ntn_repulsions(const BaseSysData& baseData, const Parameters& pars, bool Hessian, const long& timeStep);


    void nts_find_closest_approach(const int& node, const int& segment,const int& masterMesh, const BaseSysData& baseData, const Parameters& pars);

    void affine_axial_shearing(const BaseSysData& baseData, const Parameters& pars, double strain);
    void special_localized_deformation(const BaseSysData& baseData, const Parameters& pars,const double& gammaX, const double& gammaY, const std::vector<int>& targetNodes);
    void affine_compression(const BaseSysData& baseData, const Parameters& pars, double strain);
    void hold(const BaseSysData& baseData, const Parameters& pars);
    void dump_global_data(const Parameters& pars, const long& timeStep, std::string mode, std::string purpose);  //mode: "w" for writing or "a" for appending. purpose: "i" for inspection or "f" for final results
    

    
    // Variables and functions needed for the Hessian
    
    Eigen::VectorXd WmPrimePrimeJSquared; // Wm is mixing energy. PrimePrime is the second deriviative w/r/t areaRatio (aka J)
    Eigen::VectorXd prefactorA; // = (NkT + Wm" J^2)
    Eigen::VectorXd prefactorB;  // = (NkT - Wm'J)
    Eigen::VectorXd prefactorC; // = (Wm" J^2 + Wm' J)
    Eigen::VectorXd Kxxxx;
    Eigen::VectorXd Kxxxy;
    Eigen::VectorXd Kxxyx;
    Eigen::VectorXd Kxxyy;
    Eigen::VectorXd Kxyxy;
    Eigen::VectorXd Kxyyx;
    Eigen::VectorXd Kxyyy;
    Eigen::VectorXd Kyxyx;
    Eigen::VectorXd Kyxyy;
    Eigen::VectorXd Kyyyy;
    Eigen::SparseMatrix<double> KMxxjx;
    Eigen::SparseMatrix<double> KMxxjy;
    Eigen::SparseMatrix<double> KMxyjx;
    Eigen::SparseMatrix<double> KMxyjy;
    Eigen::SparseMatrix<double> KMyxjx;
    Eigen::SparseMatrix<double> KMyxjy;
    Eigen::SparseMatrix<double> KMyyjx;
    Eigen::SparseMatrix<double> KMyyjy;
    Eigen::SparseMatrix<double> Hixjx;
    Eigen::SparseMatrix<double> Hixjy;
    Eigen::SparseMatrix<double> Hiyjx;
    Eigen::SparseMatrix<double> Hiyjy;
    Eigen::SparseMatrix<double> InvHixjx;
    Eigen::SparseMatrix<double> InvHixjy;
    Eigen::SparseMatrix<double> InvHiyjx;
    Eigen::SparseMatrix<double> InvHiyjy;
    Eigen::VectorXd affineForceX;
    Eigen::VectorXd affineForceY;
    Eigen::VectorXd nonAffineVX;
    Eigen::VectorXd nonAffineVY;
    Eigen::SparseMatrix<double> Hessian;
    Eigen::VectorXd augmentedAffineF;
    Eigen::VectorXd augmentedNonAffineV;
//    Eigen::SparseMatrix<double> contactHixjx;
//    Eigen::SparseMatrix<double> contactHixjy;
//    Eigen::SparseMatrix<double> contactHiyjx;
//    Eigen::SparseMatrix<double> contactHiyjy;
    
//    Eigen::VectorXd FXref;
//    Eigen::VectorXd FYref;
//    Eigen::VectorXd DeltaForceXRatio;
//    Eigen::VectorXd DeltaForceYRatio;
//    Eigen::VectorXd refX;
//    Eigen::VectorXd refY;
//    Eigen::VectorXd HessianFX;
//    Eigen::VectorXd HessianFY;
    
    void calculate_monolithic_stiffness_tensor(const Parameters& pars);
    void calculate_the_hessian(const Parameters& pars);
    void add_d1_contributions_to_the_hessian(double penaltyStifness, double xs,double ys,double x0,double y0,double x1,double y1, int snode, int node0, int node1, const BaseSysData& baseData);
    void add_d2_contributions_to_the_hessian(double penaltyStifness, double xs,double ys,double xm,double ym, int snode, int mnode, const BaseSysData& baseData);
    void fill_augmented_Hessian();
    

};

#endif /* Configuration_hpp */
