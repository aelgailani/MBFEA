
#ifndef Parameters_hpp
#define Parameters_hpp

#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
class Parameters {
private:
    template <typename T>
    void printit(const std::string name,const T v) const;

    //    template <typename T>
//    void print2(const std::string name,const T v, std::ofstream myfile) const

public:
    Parameters(std::string& inputFileName, std::string& inputRestartFolder, std::string& runMode, std::string& restartStepOverwrite);
    double kTOverOmega;
    double NkT;
    double chi;
    double verletCellCutoff;
    double initialStretch;
    long dumpEvery;
    long splitDataEvery;
    long startingTimeStep;
    std::string outputFolderName;
    std::string trianglesFileName;
    std::string surfaceNodesFileName;
    std::string initialNodesFileName;
    std::string boundaryType;
    double imagesMargin;
    double initTopPos;
    double initBotPos;
    double initRightPos;
    double initLeftPos;
    double deformationRate;
    double targetPhi;
    double dt;
    std::string wallStyle;
    std::string solver;
    double HWallStiffness;
    double PLWallEnergyScale;
    double PLWallLJScale;
    double ntsHarmonicPenaltyStiffness;
    double Ap;
    std::string runMode;
    std::string startingMode;
    std::string restartFile;
    double targetShear;
    double maxForceTol;
    double FIRE_dtmax;
    long FIRE_N_positive_min;
    double FIRE_finc;
    double FIRE_fdec;
    double FIRE_alpha_start;
    double FIRE_falpha;
    double FIRE_dt_start;
    double FIRE_RTolerance;
    long FIRE_N_negative_max;
    double FIRE_dtmin;
    bool FIRE_intialdelay;
    long FIRE_Nmax;
    int FIRE_numStrainSteps;
    long FIRE_startingStrainStep;
    double  gammaX;
    double  gammaY;
    int segmentCellMethod;
    std::vector<int> targetNodes;
    bool dumpPeriodicImagesXY;
    bool calculateHessian;
    bool callPythonPlot;
    bool identifyAndDumbFacets;
    bool reversibleMasterSlaveRole;
    std::string contactMethod;
    std::string ntsPenaltyMethod;
    double ntnRepulseEnergy;
    double ntsPowerlawRepulseEnergy;
    double ntsPowerlawLjScale;
    double ntnLjScale;
    double ntnRcutoff;
    int integrator;
    int gntn_NGhostNodes;
    void print_to_console(void) const;
    
};

#endif /* Parameters_hpp */
