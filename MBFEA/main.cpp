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
    std::string sartingMode = "none";
    std::string runMode = "new";
    std::string restartStep = "none" ;  //restarting step
    
    // This section is for passing parameters as arguments. NOT finishid yet so just use the parameters input file for now
    if (argc>1) {
        for (int i=1;i<argc; i+=2){
            if (std::string(argv[i]) == "--parsfile" || std::string(argv[i])== "-pf") {
                inputFileName = std::string(argv[i+1]);
                std::cout << "Input file read is "<<inputFileName<< std::endl;
            }
        }
    }
//            else if (std::string(argv[i]) == "--compressrestartstep" || std::string(argv[i])== "-cr") {
//                restartStep = std::string(argv[i+1]);
//                std::cout << "restarting step is  "<<restartStep<< std::endl;
//                sartingMode = "restart";
//                std::string runMode = "compression";
//            }else if (std::string(argv[i]) == "--contshearrestartstep" || std::string(argv[i])== "-csr") {
//                restartStep = std::string(argv[i+1]);
//                std::cout << "restarting step is  "<<restartStep<< std::endl;
//                sartingMode = "restart";
//                runMode = "continuousShear";
//            }else if (std::string(argv[i]) == "--stepshearrestartstep" || std::string(argv[i])== "-ssr") {
//                restartStep = std::string(argv[i+1]);
//                std::cout << "restarting step is  "<<restartStep<< std::endl;
//                sartingMode = "restart";
//                runMode = "stepShear";
//            }else {
//                std::cout << "unexpected argument buddy ! "<<std::string(argv[i])<< std::endl;;
//                exit(1);}
//        }
//    }

    
    //Read parameters
    const Parameters pars(inputFileName, runMode, restartStep);
    pars.print_to_console();
    
    
    //Open/create directory to dump the outputs
    if (pars.startingMode == "new" || (pars.startingMode == "restart" && pars.runMode == "compress") )
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
    }else if (pars.startingMode == "restart" && pars.runMode == "stepShear"){
        DIR* dir1 = opendir(pars.outputFolderName.c_str());
        DIR* dir2 = opendir((pars.outputFolderName +"/step-"+std::to_string(pars.startingTimeStep)).c_str());
        if (dir1)
        {
            closedir(dir1);
        } else if (ENOENT == errno){
            const int dir_err = mkdir(pars.outputFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if (-1 == dir_err) {
                printf("Error creating shearing directory!");
                exit(1);
            }
        }
    
        if (dir2){
            closedir(dir2);
        } else if (ENOENT == errno){
            const int dir_err = mkdir((pars.outputFolderName +"/step-"+std::to_string(pars.startingTimeStep)).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if (-1 == dir_err) {
                printf("Error creating step directory!");
                exit(1);
            }
        }
        }else if (pars.startingMode == "restart" && pars.runMode == "continuousShear"){
            DIR* dir1 = opendir(pars.outputFolderName.c_str());
            if (dir1)
            {
                closedir(dir1);
            } else if (ENOENT == errno){
                const int dir_err = mkdir(pars.outputFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                if (-1 == dir_err) {
                    printf("Error creating shearing directory!");
                    exit(1);
                }
            }
        
            
            
        }else{printf("Please check your runMode and startingMode in the pars file. This is 'spelling-sensitive'!");
            exit(1);
            
        }
        
        
        
    

    //Fill in input data
    BaseSysData baseData(pars);
    baseData.dump_augmented_surface_meshes(pars);   // uncomment it for debugging purposes if you like to use it
    //Create system confgiuration
    Configuration mainSys(baseData, pars);
    
   
    
    int timeStep = pars.startingTimeStep;
    int stage = 0; // a dummy varaiable to be used in shearing mode
    mainSys.dump_global_data(pars, 'w', 'i');

    
    ///////////////  Main loop
    
    if (pars.solver=="GD") {
        GD_solver(baseData,pars,timeStep,stage, mainSys);
    }else if (pars.solver=="FIRE"){
        FIRE_solver(baseData,pars,timeStep,stage, mainSys);
    }
    
    
    return 0;
    
}

