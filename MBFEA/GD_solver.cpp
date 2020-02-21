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

void GD_solver(const BaseSysData& baseData, const Parameters& pars, int timeStep, int stage, Configuration& mainSys){
//    auto t0 = std::chrono::high_resolution_clock::now();
    if (pars.runMode=="compress"){  //  compressing ***********************************************************************************************
        while (1)
        {
            
            //for debugging only (done its purpose but could be needed)
//            if (timeStep>500) {
//
//                mainSys.compute_forces_PBC(baseData, pars, timeStep, 1, 0);
//
//                mainSys.FXref = mainSys.forceX;
//                mainSys.FYref = mainSys.forceY;
//                mainSys.refY = mainSys.curPosY;
//                mainSys.refX = mainSys.curPosX;
//
//                mainSys.curPosX(0)= mainSys.curPosX(0)+0.000001;
//                mainSys.curPosY(18)= mainSys.curPosY(18)+0.000001;
//                mainSys.curPosX(5)= mainSys.curPosX(5)+0.000001;
//                mainSys.curPosY(120)= mainSys.curPosY(120)+0.000001;
//                mainSys.curPosY(4)= mainSys.curPosY(4)+0.0003;
//                mainSys.curPosX(5)= mainSys.curPosX(5)+0.0001;
//                mainSys.curPosY(5)= mainSys.curPosY(5)+0.0005;
//                mainSys.curPosX(70)= mainSys.curPosX(70)+0.0001;
//                mainSys.curPosY(70)= mainSys.curPosY(70)+0.0005;
//
//                mainSys.HessianFX = mainSys.Hixjx * (mainSys.curPosX - mainSys.refX) + mainSys.Hixjy * (mainSys.curPosY - mainSys.refY);
//                mainSys.HessianFY = mainSys.Hiyjx * (mainSys.curPosX - mainSys.refX) + mainSys.Hiyjy * (mainSys.curPosY - mainSys.refY);
//
//                mainSys.compute_forces_PBC(baseData, pars, timeStep, 1, 0);
//
//                mainSys.DeltaForceXRatio = (mainSys.forceX - mainSys.FXref).array()/mainSys.HessianFX.array();
//                mainSys.DeltaForceYRatio = (mainSys.forceY - mainSys.FYref).array()/mainSys.HessianFY.array();
//
//                std::cout << "DeltaForceXRatio \n " <<mainSys.DeltaForceXRatio<< std::endl;
//                std::cout << "DeltaForceYRatio \n " <<mainSys.DeltaForceYRatio<< std::endl;
//
//                mainSys.dump_per_node(baseData, pars, timeStep);
//                mainSys.dump_per_ele(baseData, pars,timeStep);
//
//                mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
//
//
//                plotWithPython(timeStep);
//
//                exit(1);
//            }
            
//            auto t1 = std::chrono::high_resolution_clock::now();
            std::cout << timeStep << std::endl;
            
            if (timeStep * pars.deformationRate * pars.dt <= pars.maxCompression)
            {
                mainSys.compress(baseData, pars, pars.deformationRate * pars.dt);
            }

            // Take an Euler step
            if (pars.boundaryType == "walls"){
                mainSys.compute_forces_PBC(baseData, pars, timeStep, 1, 1);
            }else if (pars.boundaryType == "periodic"){
                mainSys.compute_forces_PBC(baseData, pars, timeStep, 1, 1);
            }
            
            mainSys.curPosX = mainSys.curPosX.array() + mainSys.forceX.array() * pars.dt;
            mainSys.curPosY = mainSys.curPosY.array() + mainSys.forceY.array() * pars.dt;
            mainSys.displacementSinceLastStep = ((mainSys.curPosX.array() - mainSys.curPosXAtLastStep.array()).pow(2)+(mainSys.curPosY.array()-mainSys.curPosYAtLastStep.array()).pow(2)).pow(0.5);
            if (mainSys.displacementSinceLastStep.maxCoeff() >= pars.verletCellCutoff){
                // updated curPos AtLastGridUpdate
                mainSys.curPosXAtLastStep = mainSys.curPosX;
                mainSys.curPosYAtLastStep = mainSys.curPosY;
            }

            // Postporcesseing calculations
            mainSys.update_post_processing_data(baseData, pars);
            
            
            //dump
            mainSys.dump_global_data(pars, 'a', 'i');
            if (timeStep % pars.dumpEvery == 0) {
//                mainSys.check_force_energy_consistency(baseData, pars);
                mainSys.dump_per_node(baseData, pars, timeStep);
                mainSys.dump_per_ele(baseData, pars,timeStep);
                if (pars.dumpPeriodicImagesXY){
                    mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
                }
                if (pars.callPythonPlot) {
                    plotWithPython(timeStep);
                }
            }

//            auto t2 = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double> elapsed = t2 - t1;
//            std::chrono::duration<double> elapsed0 = t2 - t0;
            
            timeStep++;
            std::cout << "maxForce  " << mainSys.maxR << std::endl;
            std::cout << "meanForce  " << mainSys.avgR << std::endl;
            std::cout << "maxDisplacement  " << mainSys.displacementSinceLastStep.maxCoeff()<< std::endl;
            std::cout << "phi  " << mainSys.phi << std::endl;
            std::cout << "deltaEnergy  " << mainSys.deltaTotEnergy << std::endl;
            std::cout << "deltaEnergy/dt  " << mainSys.deltaTotEnergy/pars.dt << std::endl;
            std::cout << "L2NormResidual  " << mainSys.L2NormResidual << std::endl;
//            std::cout << "elapsed time per step:  " << elapsed.count() << std::endl;
//            std::cout << "elapsed total:  " << elapsed0.count() << std::endl;
            std::cout << "\n" << std::endl;
            
            if ( (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) > 1E10  || (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) < 1E-10 || mainSys.maxR>50000.0){
                std::cout << "Foce condition met !" << std::endl;
                break;
            }
            if ( isnan(mainSys.areaRatio.sum()) || isnan(mainSys.forceX.sum()) ||  isnan(mainSys.forceY.minCoeff()) ||  isnan(mainSys.maxR) ){
                std::cout << "System blew up !" << std::endl;
                break;
            }
        }
    }else if (pars.runMode=="stepShear"){  //  shearing ***********************************************************************************************
        mainSys.maxR = 9999;
        while (1)
        {
           
//            auto t1 = std::chrono::high_resolution_clock::now();
            std::cout << timeStep << std::endl;
            
            
            
            if (timeStep * pars.deformationRate * pars.dt < pars.maxShear )
            {
                mainSys.shear(baseData, pars, pars.deformationRate * pars.dt);
                
            }else if (mainSys.maxR >= pars.maxForceTol)
            {
                mainSys.hold(baseData, pars);
                
            }else if (stage==0){
                
                mainSys.shear(baseData, pars, pars.maxShear);
                stage++;
                
            }else if (stage==1){
                stage++;
            }
            
            
            // Take an Euler step
            if (pars.boundaryType == "walls"){
                mainSys.compute_forces_PBC(baseData, pars, timeStep, 1, 1);
            }else if (pars.boundaryType == "periodic"){
                mainSys.compute_forces_PBC(baseData, pars, timeStep, 1, 1);
            }
            
            mainSys.curPosX = mainSys.curPosX.array() + mainSys.forceX.array() * pars.dt;
            mainSys.curPosY = mainSys.curPosY.array() + mainSys.forceY.array() * pars.dt;
            mainSys.displacementSinceLastStep = ((mainSys.curPosX.array() - mainSys.curPosXAtLastStep.array()).pow(2)+(mainSys.curPosY.array()-mainSys.curPosYAtLastStep.array()).pow(2)).pow(0.5);
            if (mainSys.displacementSinceLastStep.maxCoeff() >= pars.verletCellCutoff){
                // updated curPos AtLastGridUpdate
                mainSys.curPosXAtLastStep = mainSys.curPosX;
                mainSys.curPosYAtLastStep = mainSys.curPosY;
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
            
//            auto t2 = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double> elapsed = t2 - t1;
//            std::chrono::duration<double> elapsed0 = t2 - t0;
            timeStep++;
            std::cout << "stage  " << stage << std::endl;
            std::cout << "maxForce  " << mainSys.maxR << std::endl;
            std::cout << "meanForce  " << mainSys.avgR << std::endl;
            std::cout << "maxDisplacement  " << mainSys.displacementSinceLastStep.maxCoeff() << std::endl;
//            std::cout << "elapsed time per step:  " << elapsed.count() << std::endl;
//            std::cout << "elapsed total:  " << elapsed0.count() << std::endl;
            std::cout << "\n" << std::endl;
            
            if (stage==2){
                std::cout << "Done !" << std::endl;
                break;
            }
            if ( (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) > 1E10  || (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) < 1E-10 || mainSys.maxR>50000.0){
                std::cout << "Foce condition met !" << std::endl;
                break;
            }
            if ( isnan(mainSys.areaRatio.sum()) || isnan(mainSys.forceX.sum()) ||  isnan(mainSys.forceY.minCoeff()) ||  isnan(mainSys.maxR) ){
                std::cout << "System blew up !" << std::endl;
                break;
            }
        }
    }else if (pars.runMode=="special"){  //  compressing ***********************************************************************************************
        while (1)
        {
            
            std::cout << timeStep << std::endl;
        
            if (timeStep==0){
                mainSys.special_localized_deformation(baseData, pars, pars.gammaX, pars.gammaY, pars.targetNodes);
            }
         
            
            // Take an Euler step
            if (pars.boundaryType == "walls"){
                mainSys.compute_forces_PBC(baseData, pars, timeStep, 0, 1);
            }else if (pars.boundaryType == "periodic"){
                mainSys.compute_forces_PBC(baseData, pars, timeStep, 0, 1);
            }
            for(int nodeID: pars.targetNodes){
                mainSys.forceX(nodeID) = 0;
                mainSys.forceY(nodeID) = 0;
            }
            mainSys.curPosX = mainSys.curPosX.array() + mainSys.forceX.array() * pars.dt;
            mainSys.curPosY = mainSys.curPosY.array() + mainSys.forceY.array() * pars.dt;
//            mainSys.displacementSinceLastStep = ((mainSys.curPosX.array() -               mainSys.curPosXAtLastStep.array()).pow(2)+(mainSys.curPosY.array()-                 mainSys.curPosYAtLastStep.array()).pow(2)).pow(0.5);
//            if (mainSys.displacementSinceLastStep.maxCoeff() >= pars.verletCellCutoff){
//                // updated curPos AtLastGridUpdate
//                mainSys.curPosXAtLastStep = mainSys.curPosX;
//                mainSys.curPosYAtLastStep = mainSys.curPosY;
//            }
//

            
            // Postporcesseing calculations
            mainSys.update_post_processing_data(baseData, pars);
            
            
            //dump
            mainSys.dump_global_data(pars, 'a', 'i');
//            if (timeStep % pars.dumpEvery == 0) {
//                mainSys.check_force_energy_consistency(baseData, pars);
//                mainSys.dump_per_node(baseData, pars, timeStep);
//                mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
//                mainSys.dump_per_ele(baseData, pars,timeStep);
//                plotWithPython(timeStep);
//            }

            
            timeStep++;
//            std::cout << "maxForce  " << mainSys.maxR << std::endl;
//            std::cout << "maxDisplacement  " << mainSys.displacementSinceLastStep.maxCoeff()<< std::endl;
//            std::cout << "phi  " << mainSys.phi << std::endl;
//            std::cout << "deltaEnergy  " << mainSys.deltaTotEnergy << std::endl;
//            std::cout << "deltaEnergy/dt  " << mainSys.deltaTotEnergy/pars.dt << std::endl;
//            std::cout << "L2NormResidual  " << mainSys.L2NormResidual << std::endl;

            std::cout << "\n" << std::endl;
            
            if ( (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) > 1E10  || (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) < 1E-10 || mainSys.maxR>50000.0){
                std::cout << "Foce condition met !" << std::endl;
                break;
            }
            if ( isnan(mainSys.areaRatio.sum()) || isnan(mainSys.forceX.sum()) ||  isnan(mainSys.forceY.minCoeff()) ||  isnan(mainSys.maxR) ){
                std::cout << "System blew up !" << std::endl;
                break;
            }
        }
    }else if (pars.runMode=="continuousShear"){  //  shearing ***********************************************************************************************
        mainSys.maxR = 9999;
        while (1)
        {
        
//            auto t1 = std::chrono::high_resolution_clock::now();
            std::cout << timeStep << std::endl;
            
            
            
            if (timeStep * pars.deformationRate * pars.dt < pars.maxShear )
            {
                mainSys.shear(baseData, pars, pars.deformationRate * pars.dt);
               
            }else if (timeStep * pars.deformationRate * pars.dt >= pars.maxShear )
            {
                mainSys.hold(baseData, pars);
                
            }
            
            // Take an Euler step
//            if (pars.boundaryType == "walls"){
//                mainSys.compute_forces_PBC(baseData, pars, timeStep, 1, 1);
//            }else
            if (pars.boundaryType == "periodic"){
                mainSys.compute_forces_PBC(baseData, pars, timeStep, 1, 1);
            }
            
            mainSys.curPosX = mainSys.curPosX.array() + mainSys.forceX.array() * pars.dt;
            mainSys.curPosY = mainSys.curPosY.array() + mainSys.forceY.array() * pars.dt;
            mainSys.displacementSinceLastStep = ((mainSys.curPosX.array() - mainSys.curPosXAtLastStep.array()).pow(2)+(mainSys.curPosY.array()-mainSys.curPosYAtLastStep.array()).pow(2)).pow(0.5);
            if (mainSys.displacementSinceLastStep.maxCoeff() >= pars.verletCellCutoff){
                // updated curPos AtLastGridUpdate
                mainSys.curPosXAtLastStep = mainSys.curPosX;
                mainSys.curPosYAtLastStep = mainSys.curPosY;
            }
            
            // Postporcesseing calculations
            mainSys.update_post_processing_data(baseData, pars);
            
            
            // Dump data
            mainSys.dump_global_data(pars, 'a', 'i');
            
            
            if (timeStep % pars.dumpEvery == 0) {
                mainSys.dump_per_node(baseData, pars, timeStep);
                mainSys.dump_per_ele(baseData, pars,timeStep);
                if (pars.dumpPeriodicImagesXY){
                    mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
                }
                if (pars.callPythonPlot) {
                    plotWithPython(timeStep);
                }
                    
            }
            
//            auto t2 = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double> elapsed = t2 - t1;
//            std::chrono::duration<double> elapsed0 = t2 - t0;
            timeStep++;
            
            std::cout << "maxForce  " << mainSys.maxR << std::endl;
            std::cout << "meanForce  " << mainSys.avgR << std::endl;
//            std::cout << "maxDisplacement  " << mainSys.displacementSinceLastStep.maxCoeff() << std::endl;
//            std::cout << "elapsed time per step:  " << elapsed.count() << std::endl;
//            std::cout << "elapsed total:  " << elapsed0.count() << std::endl;
            std::cout << "\n" << std::endl;
            
           if (timeStep * pars.deformationRate * pars.dt >= pars.maxShear )
           {
               std::cout << "Done successfully!" << std::endl;
               exit(1);
           }
            if ( (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) > 1E10  || (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) < 1E-10 || mainSys.maxR>50000.0){
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
