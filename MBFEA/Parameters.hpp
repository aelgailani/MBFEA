
#ifndef Parameters_hpp
#define Parameters_hpp

#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
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
    double HWallStiffness;
    double PLWallEnergyScale;
    double PLWallLJScale;
    double penaltyStiffness;
    double Ap;
    std::string runMode;
    std::string startingMode;
    std::string restartFile;
    double shearStep;
    double maxForceTol;
    void print_to_console(void) const;
    
    
};

#endif /* Parameters_hpp */
