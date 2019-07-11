#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unordered_set>
#include <cmath>
#include <sys/stat.h>
#include <dirent.h>
#include <iomanip>
#include <chrono>
#include "Parameters.hpp"
#include "visualization.hpp"
#include "BaseSysData.hpp"
#include "Configuration.hpp"
#include "solvers.hpp"

int main(int argc, char* argv[])
{
   
    //Assign input file name to default if not passed
    std::string inputFileName = "setParameters.txt";  //default file name
    std::string sartingMode = "new";
    std::string runMode = "compress";
    std::string restartStep = "none" ;  //restarting step for shearing
    if (argc>1) {
        for (int i=1;i<argc; i+=2){
            if (std::string(argv[i]) == "--parsfile" || std::string(argv[i])== "-pf") {
                inputFileName = std::string(argv[i+1]);
                std::cout << "Input file read is "<<inputFileName<< std::endl;
            }else if (std::string(argv[i]) == "--compressrestartstep" || std::string(argv[i])== "-cr") {
                restartStep = std::string(argv[i+1]);
                std::cout << "restarting step is  "<<restartStep<< std::endl;
                sartingMode = "restart";
            }else if (std::string(argv[i]) == "--shearrestartstep" || std::string(argv[i])== "-sr") {
                restartStep = std::string(argv[i+1]);
                std::cout << "restarting step is  "<<restartStep<< std::endl;
                sartingMode = "restart";
                runMode = "shear";
            }else {
                std::cout << "unexpected argument buddy ! "<<std::string(argv[i])<< std::endl;;
                exit(1);}
        }
    }

    
    //Read parameters
    const Parameters pars(inputFileName, runMode, restartStep);
    pars.print_to_console();
    
    
    //Open/create directory to dump the outputs
    if (sartingMode == "new")
    {    DIR* dir = opendir(pars.outputFolderName.c_str());
        if (dir)
        {
            closedir(dir);
        } else if (ENOENT == errno){
            const int dir_err = mkdir(pars.outputFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if (-1 == dir_err) {
                printf("Error creating main run directory!");
                exit(1);
            }
        }
    }else if (sartingMode == "restart" && runMode == "shear"){
        DIR* dir1 = opendir("shearing");
        DIR* dir2 = opendir(pars.outputFolderName.c_str());
        if (dir1)
        {
            closedir(dir1);
        } else if (ENOENT == errno){
            const int dir_err = mkdir("shearing", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if (-1 == dir_err) {
                printf("Error creating shearing directory!");
                exit(1);
            }
        }
    
        if (dir2){
            closedir(dir2);
        } else if (ENOENT == errno){
            const int dir_err = mkdir(pars.outputFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if (-1 == dir_err) {
                printf("Error creating step directory!");
                exit(1);
            }
        }
        }
        
        
        
    

    //Fill in input data
    BaseSysData baseData(pars);

    //Create system confgiuration
    Configuration mainSys(baseData, pars);
    
    
    int timeStep = pars.startingTimeStep;
    int stage = 0; // a dummy varaiable to be used in shearing mode
    mainSys.dump_global_data(pars, 'w', 'i');

    ///////////////  Main loop
    
    if (pars.solver=="GD") {
        GD_solver(baseData,pars,timeStep,stage, mainSys);
    }
    
    
    return 0;
    
}

