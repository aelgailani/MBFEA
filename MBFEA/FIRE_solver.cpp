//
//  GD_solver.cpp
//  MBFEA
//
//  Created by Ahmed Elgailani on 7/10/19.
//  Copyright © 2019 Ahmed Elgailani. All rights reserved.
//
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
#include "solvers.hpp"
#include "Configuration.hpp"
#include "visualization.hpp"

void FIRE_solver(const BaseSysData& baseData, const Parameters& pars, int timeStep, int stage, Configuration& mainSys){
    auto t0 = std::chrono::high_resolution_clock::now();

    double FIRE_alpha =  pars.FIRE_alpha_start;
    double FIRE_N = 0;
    double FIRE_dt = pars.FIRE_dt_start;
//    double RTolerance;
    if (pars.runMode=="compress"){  //  compressing ***********************************************************************************************
        for (int strainStep = 0; strainStep<= pars.numStrainSteps ; strainStep++) {
            
            auto t1 = std::chrono::high_resolution_clock::now();
            
            
            mainSys.compress(baseData, pars, pars.maxCompression * float(strainStep)/float(pars.numStrainSteps));

            while (1)
            {
                std::cout << "strainStep  " << strainStep<<"   out of  " << pars.numStrainSteps <<std::endl;
                std::cout << "energy  " << mainSys.totalEnergy <<std::endl;
                std::cout << timeStep << std::endl;
                // Take an Euler step
                if (pars.boundaryType == "walls"){
                    mainSys.compute_forces_PBC(baseData, pars, timeStep);
                }else if (pars.boundaryType == "periodic"){
                    mainSys.compute_forces_PBC(baseData, pars, timeStep);
                }
                
                double power = mainSys.forceX.dot(mainSys.velocityX)+ mainSys.forceY.dot(mainSys.velocityY);
                std::cout << "power   " << power <<std::endl;
                
                if (power <=0) {
                    FIRE_dt *= pars.FIRE_fdec;
                    mainSys.velocityX.fill(0);
                    mainSys.velocityY.fill(0);
                    FIRE_alpha = pars.FIRE_alpha_start;
                    FIRE_N = timeStep;
                }else if ( (timeStep - FIRE_N) > pars.FIRE_Nmin){
                    FIRE_dt = fmin(FIRE_dt * pars.FIRE_finc, pars.FIRE_dtmax);
                    FIRE_alpha *= pars.FIRE_falpha;
                    Eigen::VectorXd force_magnitude = (mainSys.forceX.array().pow(2)+mainSys.forceY.array().pow(2)).pow(0.5);
                    Eigen::VectorXd velocity_magnitude = (mainSys.velocityX.array().pow(2)+mainSys.velocityY.array().pow(2)).pow(0.5);
                    mainSys.velocityX = (1 - FIRE_alpha)*mainSys.velocityX.array()+FIRE_alpha*velocity_magnitude.array()*mainSys.forceX.array()/force_magnitude.array();
                    mainSys.velocityY = (1 - FIRE_alpha)*mainSys.velocityY.array()+FIRE_alpha*velocity_magnitude.array()*mainSys.forceY.array()/force_magnitude.array();
                }
                std::cout << "FIRE_dt   " << FIRE_dt <<std::endl;
                mainSys.curPosX += mainSys.velocityX*FIRE_dt;
                mainSys.curPosY += mainSys.velocityY*FIRE_dt;
                mainSys.velocityX += mainSys.forceX*FIRE_dt;
                mainSys.velocityY += mainSys.forceY*FIRE_dt;


                mainSys.displacementSinceLastGridUpdate = ((mainSys.curPosX.array() - mainSys.curPosXAtLastGridUpdate.array()).pow(2)+(mainSys.curPosY.array()-mainSys.curPosYAtLastGridUpdate.array()).pow(2)).pow(0.5);
                
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
                std::chrono::duration<double> elapsed = t2 - t1;
                std::chrono::duration<double> elapsed0 = t2 - t0;
                
                timeStep++;
                std::cout << "maxForce  " << mainSys.maxR << std::endl;
                std::cout << "maxDisplacement  " << mainSys.displacementSinceLastGridUpdate.maxCoeff()<< std::endl;
                std::cout << "phi  " << mainSys.phi << std::endl;
                std::cout << "elapsed time per step:  " << elapsed.count() << std::endl;
                std::cout << "elapsed total:  " << elapsed0.count() << std::endl;
                std::cout << "\n" << std::endl;
                
                if ( (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) > 1E10  || (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) < pars.RTolerance){
                    std::cout << "Foce condition met !" << std::endl;
                    break;
                }

            }
            
            FIRE_alpha = pars.FIRE_alpha_start;
            FIRE_N = timeStep;
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
                mainSys.compute_forces_PBC(baseData, pars, timeStep);
            }else if (pars.boundaryType == "periodic"){
                mainSys.compute_forces_PBC(baseData, pars, timeStep);
            }
            
            mainSys.curPosX = mainSys.curPosX.array() + mainSys.forceX.array() * pars.dt;
            mainSys.curPosY = mainSys.curPosY.array() + mainSys.forceY.array() * pars.dt;
            mainSys.displacementSinceLastGridUpdate = ((mainSys.curPosX.array() - mainSys.curPosXAtLastGridUpdate.array()).pow(2)+(mainSys.curPosY.array()-mainSys.curPosYAtLastGridUpdate.array()).pow(2)).pow(0.5);
            if (mainSys.displacementSinceLastGridUpdate.maxCoeff() >= pars.verletCellCutoff){
                // updated curPos AtLastGridUpdate
                mainSys.curPosXAtLastGridUpdate = mainSys.curPosX;
                mainSys.curPosYAtLastGridUpdate = mainSys.curPosY;
            }
            
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
            std::chrono::duration<double> elapsed = t2 - t1;
            std::chrono::duration<double> elapsed0 = t2 - t0;
            timeStep++;
            std::cout << "stage  " << stage << std::endl;
            std::cout << "maxForce  " << mainSys.maxR << std::endl;
            std::cout << "maxDisplacement  " << mainSys.displacementSinceLastGridUpdate.maxCoeff() << std::endl;
            std::cout << "elapsed time per step:  " << elapsed.count() << std::endl;
            std::cout << "elapsed total:  " << elapsed0.count() << std::endl;
            std::cout << "\n" << std::endl;
            
            if (stage==2){
                std::cout << "Done !" << std::endl;
                break;
            }
            if ( (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) > 1E10  || (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) < 1E-10 || mainSys.maxR>50.0){
                std::cout << "Foce condition met !" << std::endl;
                break;
            }
            if ( isnan(mainSys.areaRatio.sum()) || isnan(mainSys.forceX.sum()) ||  isnan(mainSys.forceY.minCoeff()) ||  isnan(mainSys.maxR) ){
                std::cout << "System blew up !" << std::endl;
                break;
            }
        }
    }
}

