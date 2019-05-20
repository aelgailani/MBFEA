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

    Eigen::Matrix<double,Eigen::Dynamic, 1> curPosX;
    Eigen::Matrix<double,Eigen::Dynamic, 1>  curPosY;
    Eigen::Matrix<double, Eigen::Dynamic, 1> defGradXX;
    Eigen::Matrix<double, Eigen::Dynamic, 1> defGradXY;
    Eigen::Matrix<double, Eigen::Dynamic, 1> defGradYX;
    Eigen::Matrix<double, Eigen::Dynamic, 1> defGradYY;
    Eigen::Matrix<double, Eigen::Dynamic, 1> defGradInvTransXX;
    Eigen::Matrix<double, Eigen::Dynamic, 1> defGradInvTransXY;
    Eigen::Matrix<double, Eigen::Dynamic, 1> defGradInvTransYX;
    Eigen::Matrix<double, Eigen::Dynamic, 1> defGradInvTransYY;
    Eigen::Matrix<double, Eigen::Dynamic, 1> forceX;
    Eigen::Matrix<double, Eigen::Dynamic, 1> forceY;
    Eigen::Matrix<double, Eigen::Dynamic, 1> interForceX;
    Eigen::Matrix<double, Eigen::Dynamic, 1> interForceY;
    Eigen::Matrix<double, Eigen::Dynamic, 1> PK1stressXX;
    Eigen::Matrix<double, Eigen::Dynamic, 1> PK1stressXY;
    Eigen::Matrix<double, Eigen::Dynamic, 1> PK1stressYX;
    Eigen::Matrix<double, Eigen::Dynamic, 1> PK1stressYY;
    Eigen::Matrix<double, Eigen::Dynamic, 1> CstressXX;
    Eigen::Matrix<double, Eigen::Dynamic, 1> CstressXY;
    Eigen::Matrix<double, Eigen::Dynamic, 1> CstressYX;
    Eigen::Matrix<double, Eigen::Dynamic, 1> CstressYY;
    Eigen::Matrix<double, Eigen::Dynamic, 1> fSquared;
    Eigen::Matrix<double, Eigen::Dynamic, 1> swellingPressurePerEle;
    Eigen::Matrix<double, Eigen::Dynamic, 1> elasticEnergyPerEle;
    Eigen::Matrix<double, Eigen::Dynamic, 1> mixingEnergyPerEle;
    Eigen::Matrix<double, Eigen::Dynamic, 1> internalEnergyPerEle;
    Eigen::Matrix<double, Eigen::Dynamic, 1> refArea;
    Eigen::Matrix<double, Eigen::Dynamic, 1> areaRatio;
    Eigen::SparseMatrix<double> gradX;
    Eigen::SparseMatrix<double> gradY;
    Eigen::Matrix<double, Eigen::Dynamic, 1> wallForceTop;
    Eigen::Matrix<double, Eigen::Dynamic, 1> wallForceBottom;
    Eigen::Matrix<double, Eigen::Dynamic, 1> wallForceRight;
    Eigen::Matrix<double, Eigen::Dynamic, 1> wallForceLeft;
    std::vector<int> masterSlave;
    double totalEnergy=0, internalEnergy=0, wallsEnergy=0, contactsEnergy=0, extW = 0, KWXX=0, KWYY=0,KWWV=0, KWWH=0;
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


    
    void update_post_processing_data(const BaseSysData& baseData, const Parameters& pars);
    void dump_per_node(const BaseSysData& baseData, const Parameters& pars, int& timeStep);
    void dump_per_ele(const BaseSysData& baseData, const Parameters& pars, int& timeStep);
    void compute_forces(const BaseSysData& baseData, const Parameters& pars);
    void shear(const BaseSysData& baseData, const Parameters& pars, double strain);
    void compress(const BaseSysData& baseData, const Parameters& pars, double strain);
    void hold(const BaseSysData& baseData, const Parameters& pars);
    void dump_global_data(const Parameters& pars, char mode, char purpose);  //mode: "w" for writing or "a" for appending. purpose: "i" for inspection or "f" for final results
};

#endif /* Configuration_hpp */
