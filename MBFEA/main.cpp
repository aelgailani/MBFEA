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

    //  Main loop
    auto t0 = std::chrono::high_resolution_clock::now();
    if (pars.runMode=="compress"){  //  compressing ***********************************************************************************************
        while (1)
        {
            std::cout << timeStep << std::endl;
            auto t1 = std::chrono::high_resolution_clock::now();
  

            if (timeStep * pars.deformationRate * pars.dt <= pars.maxCompression)
            {
                mainSys.compress(baseData, pars, pars.deformationRate * pars.dt);
            }

     
           
            // Take an Euler step
            if (pars.boundaryType == "walls"){
                mainSys.compute_forces_walls(baseData, pars);
            }else if (pars.boundaryType == "periodic"){
                mainSys.compute_forces_PBC(baseData, pars);
            }
            mainSys.curPosX = mainSys.curPosX.array() + mainSys.forceX.array() * pars.dt;
            mainSys.curPosY = mainSys.curPosY.array() + mainSys.forceY.array() * pars.dt;
            
            // Postporcesseing calculations
            mainSys.update_post_processing_data(baseData, pars);

         
            
            //dump
            mainSys.dump_global_data(pars, 'a', 'i');
            if (timeStep % pars.dumpEvery == 0) {
                mainSys.dump_per_node(baseData, pars, timeStep);
                mainSys.dump_per_ele(baseData, pars,timeStep);
                plotWithPython(timeStep);
            }
            auto t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> d1 = t2 - t1;
            std::chrono::duration<double> d2 = t2 - t0;
            
            timeStep++;
            std::cout << "maxForce  " << mainSys.maxR << std::endl;
            std::cout << "phi  " << mainSys.phi << std::endl;
            std::cout << "elapsed time per step:  " << d1.count() << std::endl;
            std::cout << "elapsed total:  " << d2.count() << std::endl;
            std::cout << "\n" << std::endl;
            
            if ( (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) > 1E10  || (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) < 1E-10 || mainSys.maxR>20.0){
                std::cout << "Foce condition met !" << std::endl;
                break;
            }
            if ( isnan(mainSys.areaRatio.sum()) || isnan(mainSys.forceX.sum()) ||  isnan(mainSys.forceY.minCoeff()) ||  isnan(mainSys.maxR) ){
                std::cout << "System blew up !" << std::endl;
                break;
            }
        }
    }else if (pars.runMode=="shear"){  //  shearing ***********************************************************************************************
        mainSys.maxR = 9999;
        while (1)
        {
            double initialFiniteShear = 0.0;
            auto t1 = std::chrono::high_resolution_clock::now();
            std::cout << timeStep << std::endl;

            
            
            if (timeStep * pars.deformationRate * pars.dt < initialFiniteShear )
            {
                mainSys.shear(baseData, pars, pars.deformationRate * pars.dt);
                
            }else if (mainSys.maxR >= pars.maxForceTol)
            {
                mainSys.hold(baseData, pars);
            
            }else if (stage==0){
                
                mainSys.shear(baseData, pars, pars.shearStep);
                stage++;
                
            }else if (stage==1){
                stage++;
            }
            

            


            
            // Take an Euler step
            if (pars.boundaryType == "walls"){
                mainSys.compute_forces_walls(baseData, pars);
            }else if (pars.boundaryType == "periodic"){
                mainSys.compute_forces_PBC(baseData, pars);
            }

            mainSys.curPosX = mainSys.curPosX.array() + mainSys.forceX.array() * pars.dt;
            mainSys.curPosY = mainSys.curPosY.array() + mainSys.forceY.array() * pars.dt;
            
            // Postporcesseing calculations
            mainSys.update_post_processing_data(baseData, pars);
            
            // Dump data
            mainSys.dump_global_data(pars, 'a', 'i');
            
            if (stage==0){
                
                mainSys.dump_global_data(pars, 'w' , 'f');
                mainSys.dump_global_data(pars, 'a' , 'f');
                
            }else if (stage==2){
                
                mainSys.dump_global_data(pars, 'a' , 'f');
                
            }
            
      
            if (timeStep % pars.dumpEvery == 0) {
                mainSys.dump_per_node(baseData, pars, timeStep);
                mainSys.dump_per_ele(baseData, pars,timeStep);
            }
            
            auto t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> d1 = t2 - t1;
            std::chrono::duration<double> d2 = t2 - t0;
            
            timeStep++;
            std::cout << "stage  " << stage << std::endl;
            std::cout << "maxForce  " << mainSys.maxR << std::endl;
            std::cout << "elapsed time per step:  " << d1.count() << std::endl;
            std::cout << "elapsed total:  " << d2.count() << std::endl;
            std::cout << "\n" << std::endl;
            
            if (stage==2){
                std::cout << "Done !" << std::endl;
                break;
            }
            if ( (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) > 1E10  || (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) < 1E-10 || mainSys.maxR>20.0){
                std::cout << "Foce condition met !" << std::endl;
                break;
            }
            if ( isnan(mainSys.areaRatio.sum()) || isnan(mainSys.forceX.sum()) ||  isnan(mainSys.forceY.minCoeff()) ||  isnan(mainSys.maxR) ){
                std::cout << "System blew up !" << std::endl;
                break;
            }
        }
    }
    
    
    return 0;
    
}

