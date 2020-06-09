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
#include "Parameters.hpp"
#include "standard_runs.hpp"
#include "BaseSysData.hpp"
#include "Configuration.hpp"
#include "solvers.hpp"
#include <time.h>

int main(int argc, char* argv[])
{

    //Assign input file name to default if not passed
    std::string inputFileName = "setParameters.txt";  //default file name
    std::string sartingMode = "new";
    std::string runMode = "none";
    std::string restartStep = "none" ;  //restarting step
    std::string inputRestartFolder = "none";  //default file name
    
    // This section is for passing parameters as arguments. NOT finishid yet so just use the parameters input file for now
    if (argc>1) {
        for (int i=1;i<argc; i+=2){
            if (std::string(argv[i]) == "--parsfile" || std::string(argv[i])== "-pf") {
                inputFileName = std::string(argv[i+1]);
                std::cout << "Input file read is "<<inputFileName<< std::endl;
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
            else if (std::string(argv[i]) == "--stepshearrestartstep" || std::string(argv[i])== "-ss") {
                restartStep = std::string(argv[i+1]);
                inputRestartFolder = std::string(argv[i+2]);
                std::cout << "restarting step is  "<<restartStep<< std::endl;
                sartingMode = "restart";
                runMode = "stepShear";
            }else if (std::string(argv[i]) == "--inputdirectory" || std::string(argv[i])== "-id") {
                inputRestartFolder = std::string(argv[i+1]);

            }else {
                std::cout << "unexpected argument buddy ! "<<std::string(argv[i])<< std::endl;;
                exit(1);}
        }
    
    }
    
    //Read parameters
    const Parameters pars(inputFileName, inputRestartFolder, runMode, restartStep);
    
    // print time to log file
    std::ofstream logfile;
    logfile.open (pars.runMode+"-parameters.log", std::ios_base::app);
    std::time_t t = std::time(0);   // get time now
    tm* localtm = std::localtime(&t);
    if (pars.runMode=="stepShear"){
        if (restartStep=="none") restartStep=pars.restartFile;
        logfile << std::endl << std::endl << "step  " << restartStep << "------ " << asctime(localtm) << std::endl;
    }else{
       logfile << std::endl << std::endl << "------ " << asctime(localtm) << std::endl;
    }
    
   
    logfile.close();
    pars.print_to_console();
    
    
        //Open/create directory to dump the outputs
    if (pars.startingMode == "new" || (pars.startingMode == "restart" && pars.runMode == "compress") || (pars.startingMode == "restart" && pars.runMode == "special"))
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
        DIR* dir2 = opendir((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingTimeStep))).c_str());
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
            const int dir_err = mkdir((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingTimeStep))).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
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
    
   
    
    long timeStep = pars.startingTimeStep;
    if (pars.runMode != "compress") {
        timeStep=0;
    }
   
    mainSys.dump_global_data(pars, timeStep, "write", "running");

    
    ///////////////  Main loop
    
    if (pars.runMode=="compress") {
        compress(baseData, pars, timeStep , mainSys);
    }else if (pars.runMode=="stepShear"){
        stepshear(baseData, pars, timeStep , mainSys);
    }else if (pars.runMode=="special"){
        if(pars.solver=="FIRE2"){
            shear_special_FIRE(baseData, pars, timeStep , mainSys);
        }else if (pars.solver=="GD"){
            shear_special_GD(baseData, pars, timeStep, mainSys);
        }
    }
    
    return 0;
    
}


