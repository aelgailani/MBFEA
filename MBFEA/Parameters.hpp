
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
//    void print2(const std::string name,const T v, std::ofstream myfile) const;
    
public:
    Parameters(std::string& inputFileName, std::string& runMode, std::string& restartStepOverwrite);
    double kTOverOmega;
    double NkT;
    double chi;
    double verletCellCutoff;
    double initialStretch;
    int dumpEvery;
    double startingTimeStep;
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
    double maxCompression;
    double dt;
    std::string wallStyle;
    std::string solver;
    double HWallStiffness;
    double PLWallEnergyScale;
    double PLWallLJScale;
    double penaltyStiffness;
    double Ap;
    std::string runMode;
    std::string startingMode;
    std::string restartFile;
    double maxShear;
    double maxForceTol;
    double FIRE_dtmax;
    double FIRE_Nmin;
    double FIRE_finc;
    double FIRE_fdec;
    double FIRE_alpha_start;
    double FIRE_falpha;
    double FIRE_dt_start;
    double RTolerance;
    int numStrainSteps;
    int startingStrainStep;
    double  gammaX;
    double  gammaY;
    int segmentCellMethod;
    std::vector<int> targetNodes;
    bool dumpPeriodicImagesXY;
    bool callPythonPlot;
    void print_to_console(void) const;
    

    
    
};

#endif /* Parameters_hpp */
