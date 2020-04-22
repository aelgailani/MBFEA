//
//  gd_solver.cpp
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
#include <thread>
#include "Parameters.hpp"
#include "BaseSysData.hpp"
#include "solvers.hpp"
#include "Configuration.hpp"
//#include "visualization.hpp"

void fire_solver_old(const BaseSysData& baseData, const Parameters& pars, long timeStep, int stage, Configuration& mainSys){
//    auto t0 = std::chrono::high_resolution_clock::now();
    double FIRE_alpha =  pars.FIRE_alpha_start;
    double FIRE_N = 0;
    double FIRE_dt = pars.FIRE_dt_start;
    double FIRE_prevDt = pars.FIRE_dt_start;
    
    if (pars.runMode=="compress"){  //  compressing ***********************************************************************************************
        assert(pars.startingStrainStep <= pars.numStrainSteps);
        
        for (long strainStep = pars.startingStrainStep ; strainStep<= pars.numStrainSteps ; strainStep++) {
            
            mainSys.compress(baseData, pars, pars.targetPhi * float(strainStep)/float(pars.numStrainSteps));

            while (1)
            {
                std::cout << "strainStep  " << strainStep<<"   out of  " << pars.numStrainSteps <<std::endl;
                std::cout << "energy  " << mainSys.totalEnergy <<std::endl;
                std::cout << timeStep << std::endl;
                
                // Calculate F
                if (pars.boundaryType == "walls"){
                    mainSys.compute_forces_harmonic_walls(baseData, pars, timeStep, 1, 0, pars.calculateHessian);
                }else if (pars.boundaryType == "periodic"){
                    mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, pars.calculateHessian);
                }
                
                mainSys.velocityX += mainSys.forceX*FIRE_dt;
                mainSys.velocityY += mainSys.forceY*FIRE_dt;
                
                double power = mainSys.forceX.dot(mainSys.velocityX)+ mainSys.forceY.dot(mainSys.velocityY);
                std::cout << "power   " << power <<std::endl;
                
//                if (timeStep==152716 || timeStep== 152717 || timeStep== 152715){
//                    mainSys.dump_per_node(baseData, pars, timeStep);
//                    mainSys.dump_per_ele(baseData, pars,timeStep);
//                    plotWithPython(timeStep);
//                    std::this_thread::sleep_for (std::chrono::milliseconds(700));
//                }
                
                // make sure nothing is screwed up
                
                
                if (isnan(power) || isnan(mainSys.forceX.sum()) || isnan(mainSys.forceY.sum()) || isnan(mainSys.areaRatio.sum()) || isnan(mainSys.velocityX.sum()) || isnan(mainSys.velocityY.sum()) || mainSys.areaRatio.minCoeff()<= 1.0 ){
//                    mainSys.dump_per_node(baseData, pars, timeStep);
//                    mainSys.dump_per_ele(baseData, pars,timeStep);
//                    plotWithPython(timeStep);
                    mainSys.curPosX =  mainSys.curPosXAtLastStep;
                    mainSys.curPosY = mainSys.curPosYAtLastStep;
                    mainSys.velocityX = mainSys.prevVelocityX;
                    mainSys.velocityY = mainSys.prevVelocityY;
//                    mainSys.forceX.fill(0);
//                    mainSys.forceY.fill(0);
                    FIRE_alpha = pars.FIRE_alpha_start;
                    FIRE_N = timeStep;
                    FIRE_dt = FIRE_prevDt * pars.FIRE_fdec *0.2;
                    FIRE_prevDt = FIRE_dt;
                    mainSys.curPosX += mainSys.velocityX*FIRE_dt;
                    mainSys.curPosY += mainSys.velocityY*FIRE_dt;
                    continue;

                }else
                if (power <=0) {
                    FIRE_dt *= pars.FIRE_fdec;
                    mainSys.velocityX.fill(0);
                    mainSys.velocityY.fill(0);
                    FIRE_alpha = pars.FIRE_alpha_start;
                    FIRE_N = timeStep;
                } else {
                    if ((timeStep - FIRE_N) > pars.FIRE_Nmin){
                        FIRE_prevDt = FIRE_dt;
                        FIRE_dt = fmin(FIRE_dt * pars.FIRE_finc, pars.FIRE_dtmax);
                        FIRE_alpha *= pars.FIRE_falpha;
                    }
                    Eigen::VectorXd force_magnitude = (mainSys.forceX.array().pow(2)+mainSys.forceY.array().pow(2)).pow(0.5);
                    Eigen::VectorXd velocity_magnitude = (mainSys.velocityX.array().pow(2)+mainSys.velocityY.array().pow(2)).pow(0.5);
                    mainSys.velocityX = (1 - FIRE_alpha)*mainSys.velocityX.array()+FIRE_alpha*velocity_magnitude.array()*mainSys.forceX.array()/force_magnitude.array();
                    mainSys.velocityY = (1 - FIRE_alpha)*mainSys.velocityY.array()+FIRE_alpha*velocity_magnitude.array()*mainSys.forceY.array()/force_magnitude.array();
                }
                std::cout << "FIRE_dt   " << FIRE_dt <<std::endl;
                mainSys.curPosXAtLastStep = mainSys.curPosX;
                mainSys.curPosYAtLastStep = mainSys.curPosY;
                mainSys.prevVelocityX =  mainSys.velocityX;
                mainSys.prevVelocityY =  mainSys.velocityY;
                
                
                
                mainSys.curPosX += mainSys.velocityX*FIRE_dt;
                mainSys.curPosY += mainSys.velocityY*FIRE_dt;
                
                


                mainSys.displacementSinceLastStep = ((mainSys.curPosX.array() - mainSys.curPosXAtLastStep.array()).pow(2)+(mainSys.curPosY.array()-mainSys.curPosYAtLastStep.array()).pow(2)).pow(0.5);
                
                // Postporcesseing calculations
                mainSys.update_post_processing_data(baseData, pars);
                
                
                //dump
                mainSys.dump_global_data(pars, timeStep, "append", "running");
//                if (timeStep < 188669 && timeStep > 188660) {
//                    mainSys.dump_per_node(baseData, pars, timeStep);
//                    mainSys.dump_per_ele(baseData, pars,timeStep);
//                    plotWithPython(timeStep);
//                }

//                auto t2 = std::chrono::high_resolution_clock::now();
//                std::chrono::duration<double> elapsed = t2 - t1;
//                std::chrono::duration<double> elapsed0 = t2 - t0;
                
                timeStep++;
                std::cout << "maxForce  " << mainSys.maxR << std::endl;
                std::cout << "avgForce  " << mainSys.avgR << std::endl;
//                std::cout << "maxDisplacement  " << mainSys.displacementSinceLastStep.maxCoeff()<< std::endl;
                std::cout << "phi  " << mainSys.phi << std::endl;
//                std::cout << "elapsed time per step:  " << elapsed.count() << std::endl;
//                std::cout << "elapsed total:  " << elapsed0.count() << std::endl;
                std::cout << "\n" << std::endl;
                
                if ( mainSys.maxR > 1E10  || mainSys.maxR <= pars.RTolerance){
                    std::cout << "Foce condition met !" << std::endl;
                    break;
                }

            }
//            mainSys.check_force_energy_consistency(baseData, pars);
            FIRE_alpha = pars.FIRE_alpha_start;
            FIRE_N = timeStep;
            mainSys.dump_global_data(pars, timeStep, "append", "running");
            mainSys.dump_per_node(baseData, pars, strainStep);
            mainSys.dump_per_ele(baseData, pars,strainStep);
            if (pars.dumpPeriodicImagesXY){
                mainSys.dump_per_node_periodic_images_on(baseData, pars, strainStep);
            }
           
        }
        
    }
//        else if (pars.runMode=="shear"){  //  shearing ***********************************************************************************************
//        mainSys.maxR = 9999;
//        while (1)
//        {
//            double initialFiniteShear = 0.0;
//            auto t1 = std::chrono::high_resolution_clock::now();
//            std::cout << timeStep << std::endl;
//
//
//
//            if (timeStep * pars.deformationRate * pars.dt < initialFiniteShear )
//            {
//                mainSys.shear(baseData, pars, pars.deformationRate * pars.dt);
//
//            }else if (mainSys.maxR >= pars.maxForceTol)
//            {
//                mainSys.hold(baseData, pars);
//
//            }else if (stage==0){
//
//                mainSys.shear(baseData, pars, pars.shearStep);
//                stage++;
//
//            }else if (stage==1){
//                stage++;
//            }
//
//
//            // Take an Euler step
//            if (pars.boundaryType == "walls"){
//                mainSys.compute_forces_pbc(baseData, pars, timeStep);
//            }else if (pars.boundaryType == "periodic"){
//                mainSys.compute_forces_pbc(baseData, pars, timeStep);
//            }
//
//            mainSys.curPosX = mainSys.curPosX.array() + mainSys.forceX.array() * pars.dt;
//            mainSys.curPosY = mainSys.curPosY.array() + mainSys.forceY.array() * pars.dt;
//            mainSys.displacementSinceLastStep = ((mainSys.curPosX.array() - mainSys.curPosXAtLastStep.array()).pow(2)+(mainSys.curPosY.array()-mainSys.curPosYAtLastStep.array()).pow(2)).pow(0.5);
//            if (mainSys.displacementSinceLastStep.maxCoeff() >= pars.verletCellCutoff){
//                // updated curPos AtLastGridUpdate
//                mainSys.curPosXAtLastStep = mainSys.curPosX;
//                mainSys.curPosYAtLastStep = mainSys.curPosY;
//            }
//
//            // Postporcesseing calculations
//            mainSys.update_post_processing_data(baseData, pars);
//
//
//            // Dump data
//            mainSys.dump_global_data(pars, 'a', 'i');
//
//            if (stage==0){
//
//                mainSys.dump_global_data(pars, 'w' , 'f');
//                mainSys.dump_global_data(pars, 'a' , 'f');
//
//            }else if (stage==2){
//
//                mainSys.dump_global_data(pars, 'a' , 'f');
//
//            }
//
//
//            if (timeStep % pars.dumpEvery == 0) {
//                mainSys.dump_per_node(baseData, pars, timeStep);
//                mainSys.dump_per_ele(baseData, pars,timeStep);
//            }
//
//            auto t2 = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double> elapsed = t2 - t1;
//            std::chrono::duration<double> elapsed0 = t2 - t0;
//            timeStep++;
//            std::cout << "stage  " << stage << std::endl;
//            std::cout << "maxForce  " << mainSys.maxR << std::endl;
//            std::cout << "maxDisplacement  " << mainSys.displacementSinceLastStep.maxCoeff() << std::endl;
//            std::cout << "elapsed time per step:  " << elapsed.count() << std::endl;
//            std::cout << "elapsed total:  " << elapsed0.count() << std::endl;
//            std::cout << "\n" << std::endl;
//
//            if (stage==2){
//                std::cout << "Done !" << std::endl;
//                break;
//            }
//            if ( (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) > 1E10  || (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) < 1E-10 || mainSys.maxR>50.0){
//                std::cout << "Foce condition met !" << std::endl;
//                break;
//            }
//            if ( isnan(mainSys.areaRatio.sum()) || isnan(mainSys.forceX.sum()) ||  isnan(mainSys.forceY.minCoeff()) ||  isnan(mainSys.maxR) ){
//                std::cout << "System blew up !" << std::endl;
//                break;
//            }
//        }
//    }
}

