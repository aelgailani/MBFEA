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
#include<Eigen/IterativeLinearSolvers>
#include <cmath>
#include <sys/stat.h>
#include <dirent.h>
#include <iomanip>
#include <valarray>
#include "Parameters.hpp"
#include "BaseSysData.hpp"
#include "solvers.hpp"
#include "Configuration.hpp"
//#include "visualization.hpp"

void gd_solver(const BaseSysData& baseData, const Parameters& pars, long timeStep, int stage, Configuration& mainSys){
//    auto t0 = std::chrono::high_resolution_clock::now();
    if (pars.runMode=="compress"){  //  compressing ***********************************************************************************************
        while (1)
        {
            
            //for debugging only (done its purpose but could be needed)
//            if (timeStep>1000) {
//
//                mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, 1);
//
//                mainSys.FXref = mainSys.forceX;
//                mainSys.FYref = mainSys.forceY;
//                mainSys.refY = mainSys.curPosY;
//                mainSys.refX = mainSys.curPosX;

//                mainSys.curPosX(0)= mainSys.curPosX(0)+0.000001;
//                mainSys.curPosY(18)= mainSys.curPosY(18)+0.000001;
//                mainSys.curPosX(5)= mainSys.curPosX(5)+0.000001;
//                mainSys.curPosY(120)= mainSys.curPosY(120)+0.000001;
//                mainSys.curPosY(4)= mainSys.curPosY(4)+0.00003;
//                mainSys.curPosX(10)= mainSys.curPosX(10)+0.00001;
//                mainSys.curPosY(10)= mainSys.curPosY(10)+0.00005;
//                mainSys.curPosX(40)= mainSys.curPosX(40)+0.00001;
//                mainSys.curPosY(40)= mainSys.curPosY(40)+0.00005;
//                mainSys.curPosX(70)= mainSys.curPosX(70)+0.00001;
//                mainSys.curPosY(70)= mainSys.curPosY(70)+0.00005;

//                mainSys.HessianFX = mainSys.Hixjx * (mainSys.curPosX - mainSys.refX) + mainSys.Hixjy * (mainSys.curPosY - mainSys.refY);
//                mainSys.HessianFY = mainSys.Hiyjx * (mainSys.curPosX - mainSys.refX) + mainSys.Hiyjy * (mainSys.curPosY - mainSys.refY);
//
//                mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, 1);
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
            
            
//            if (timeStep>800000){
//
//                mainSys.shear(baseData, pars, 0.001);
//                for (int i = timeStep ; i<1000000; i++){
//                mainSys.hold(baseData, pars);
//                mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, 0);
//                mainSys.curPosX = mainSys.curPosX.array() + mainSys.forceX.array() * pars.dt;
//                mainSys.curPosY = mainSys.curPosY.array() + mainSys.forceY.array() * pars.dt;
//                timeStep++;
//                std::cout << timeStep << std::endl;
//                std::cout << "maxForce  " << mainSys.maxR << std::endl;
//                std::cout << "meanForce  " << mainSys.avgR << std::endl;
//                }
//                 mainSys.dump_per_node(baseData, pars, timeStep);
//                mainSys.dump_per_ele(baseData, pars,timeStep);
//                //                                if (pars.dumpPeriodicImagesXY){
//                mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
//                //                                }
//                //                                if (pars.callPythonPlot) {
//                plotWithPython(timeStep);
//                //                                }
//                mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, 1);
//
//
//
//                mainSys.fill_augmented_Hessian();
//
//                // The following is ONLY for axial shear  case
//                mainSys.affineForceX = mainSys.Hixjx * 0.5 * (mainSys.curPosX) - mainSys.Hixjy * 0.5 *(mainSys.curPosY);
//                mainSys.affineForceY = mainSys.Hiyjx * 0.5 *(mainSys.curPosX) - mainSys.Hiyjy * 0.5 *(mainSys.curPosY);
//
//                // The following is ONLY for simple shear  case
//                mainSys.affineForceX = mainSys.Hixjx *(mainSys.curPosY);
//                mainSys.affineForceY = mainSys.Hiyjx  *(mainSys.curPosY);
//                mainSys.augmentedAffineF << mainSys.affineForceX*(-1), mainSys.affineForceY*(-1) ;
//
//
//                Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper>solver;
//                solver.setMaxIterations(1000);
//                solver.setTolerance(0.00000001);
//
//                solver.compute(mainSys.Hessian);
//                if(solver.info()!=Eigen::Success) {
//                     std::cout <<  " decomposition failed" << std::endl;
//                }
////
//                mainSys.augmentedNonAffineV = solver.solve(mainSys.augmentedAffineF);
//                std::cout << "#iterations:     " << solver.iterations() << std::endl;
//                std::cout << "estimated error: " << solver.error()      << std::endl;
//                if(solver.info()!=Eigen::Success) {
//                  std::cout <<  " solving success" << std::endl;
//                }
//
////
//                mainSys.nonAffineVX = mainSys.augmentedNonAffineV.head(baseData.numOriginalNodes);
//                mainSys.nonAffineVY = mainSys.augmentedNonAffineV.tail(baseData.numOriginalNodes);
////
//
//                mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, 1);
//                mainSys.curPosX = mainSys.curPosX.array() + mainSys.forceX.array() * pars.dt;
//                mainSys.curPosY = mainSys.curPosY.array() + mainSys.forceY.array() * pars.dt;
//                mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, 1);
//                double mu  = (mainSys.Kxyxy.array()*mainSys.refArea.array()*mainSys.areaRatio.array()).mean();
//                double mu  =((mainSys.Kxxxx.array()+mainSys.Kyyyy.array()-2*mainSys.Kxxyy.array())*mainSys.refArea.array()*mainSys.areaRatio.array()/4).mean() ;
//

//                mainSys.shear(baseData, pars, 0.001);
//               for (int i = timeStep ; i<1000000; i++){
//               mainSys.hold(baseData, pars);
//               mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, 0);
//               mainSys.curPosX = mainSys.curPosX.array() + mainSys.forceX.array() * pars.dt;
//               mainSys.curPosY = mainSys.curPosY.array() + mainSys.forceY.array() * pars.dt;
//               timeStep++;
//               std::cout << timeStep << std::endl;
//               std::cout << "maxForce  " << mainSys.maxR << std::endl;
//               std::cout << "meanForce  " << mainSys.avgR << std::endl;
//               }
//                mu -= (mainSys.affineForceX.dot(mainSys.nonAffineVX));
//                mu -= (mainSys.affineForceY.dot(mainSys.nonAffineVY));
//                mu -= (mainSys.affineForceX.dot(mainSys.InvHixjy * mainSys.affineForceY));
//                mu -= (mainSys.affineForceY.dot(mainSys.InvHiyjx * mainSys.affineForceX));
//                mu -= (mainSys.affineForceY.dot(mainSys.InvHiyjy * mainSys.affineForceY));

//                Eigen::SparseMatrix<double> yy  = (mainSys.Hixjx*mainSys.InvHixjx);
//                std::cout << "(0,0) :     " << yy.coeffRef(0, 0) << std::endl;
//                std::cout << "(0,1) :     " << yy.coeffRef(0, 1) << std::endl;
//                std::cout << "(10,140) :     " <<  (mainSys.curPosY.array()- mainSys.yMid).matrix().sum()  << std::endl;
//                std::cout << "(3,3) :     " << yy.coeffRef(3, 3) << std::endl;
//                std::cout << "Modulus :     " << mu << std::endl;
//                if (timeStep % pars.dumpEvery == 0) {
//                //                mainSys.check_force_energy_consistency(baseData, pars);
//                mainSys.dump_per_node(baseData, pars, timeStep);
//                mainSys.dump_per_ele(baseData, pars,timeStep);
//                mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
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
                mainSys.compute_forces_harmonicWalls(baseData, pars, timeStep, 1, 0, pars.calculateHessian);
            } else if (pars.boundaryType == "periodic"){
                mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, pars.calculateHessian);
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
            mainSys.dump_global_data(pars, timeStep, "append", "running");
            
            if (timeStep % pars.dumpEvery == 0) {
//                mainSys.check_force_energy_consistency(baseData, pars);
                mainSys.dump_per_node(baseData, pars, timeStep);
                mainSys.dump_per_ele(baseData, pars,timeStep);
                if (pars.dumpPeriodicImagesXY){
                    mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
                }
                if (pars.identifyAndDumbFacets) {
                    mainSys.dump_facets(baseData, pars, timeStep);
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



            if (timeStep * pars.deformationRate * pars.dt < pars.maxShear) //Here maxShear is initial nonzerovalue for shear
            {
                mainSys.shear(baseData, pars, pars.deformationRate * pars.dt);

            }else if (mainSys.maxR >= pars.maxForceTol)
            {
                mainSys.hold(baseData, pars);

            }else if (stage==0){

                mainSys.shear(baseData, pars, pars.maxShear); //here maxShear is the shear step to be taken to measure mu
                stage++;

            }else if (stage==1){
                stage++;
            }


            // Take an Euler step
            if (pars.boundaryType == "walls"){
                mainSys.compute_forces_harmonicWalls(baseData, pars, timeStep, 1, 0, pars.calculateHessian);
            }else if (pars.boundaryType == "periodic"){
                mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, pars.calculateHessian);
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
            mainSys.dump_global_data(pars, timeStep, "append", "running");

            if (stage==0){

                mainSys.dump_global_data(pars, timeStep, "write", "final");
                mainSys.dump_global_data(pars, timeStep, "append", "final");

            }else if (stage==2){

                mainSys.dump_global_data(pars, timeStep, "append", "final");

            }


            if (timeStep % pars.dumpEvery == 0) {
                mainSys.dump_per_node(baseData, pars, timeStep);
                mainSys.dump_per_ele(baseData, pars,timeStep);
                if (pars.dumpPeriodicImagesXY){
                    mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
                }
                if (pars.identifyAndDumbFacets) {
                    mainSys.dump_facets(baseData, pars, timeStep);
                }
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
    }
//    else if (pars.runMode=="special"){  //  compressing ***********************************************************************************************
//        while (1)
//        {
//
//            std::cout << timeStep << std::endl;
//
//            if (timeStep==0){
//                mainSys.special_localized_deformation(baseData, pars, pars.gammaX, pars.gammaY, pars.targetNodes);
//            }
//
//
//            // Take an Euler step
//            if (pars.boundaryType == "walls"){
//                mainSys.compute_forces_pbc(baseData, pars, timeStep, 0, 1,0);
//            }else if (pars.boundaryType == "periodic"){
//                mainSys.compute_forces_pbc(baseData, pars, timeStep, 0, 1,0);
//            }
//            for(int nodeID: pars.targetNodes){
//                mainSys.forceX(nodeID) = 0;
//                mainSys.forceY(nodeID) = 0;
//            }
//            mainSys.curPosX = mainSys.curPosX.array() + mainSys.forceX.array() * pars.dt;
//            mainSys.curPosY = mainSys.curPosY.array() + mainSys.forceY.array() * pars.dt;
////            mainSys.displacementSinceLastStep = ((mainSys.curPosX.array() -               mainSys.curPosXAtLastStep.array()).pow(2)+(mainSys.curPosY.array()-                 mainSys.curPosYAtLastStep.array()).pow(2)).pow(0.5);
////            if (mainSys.displacementSinceLastStep.maxCoeff() >= pars.verletCellCutoff){
////                // updated curPos AtLastGridUpdate
////                mainSys.curPosXAtLastStep = mainSys.curPosX;
////                mainSys.curPosYAtLastStep = mainSys.curPosY;
////            }
////
//
//
//            // Postporcesseing calculations
//            mainSys.update_post_processing_data(baseData, pars);
//
//
//            //dump
//            mainSys.dump_global_data(pars, 'a', 'i');
////            if (timeStep % pars.dumpEvery == 0) {
////                mainSys.check_force_energy_consistency(baseData, pars);
////                mainSys.dump_per_node(baseData, pars, timeStep);
////                mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
////                mainSys.dump_per_ele(baseData, pars,timeStep);
////                plotWithPython(timeStep);
////            }
//
//
//            timeStep++;
////            std::cout << "maxForce  " << mainSys.maxR << std::endl;
////            std::cout << "maxDisplacement  " << mainSys.displacementSinceLastStep.maxCoeff()<< std::endl;
////            std::cout << "phi  " << mainSys.phi << std::endl;
////            std::cout << "deltaEnergy  " << mainSys.deltaTotEnergy << std::endl;
////            std::cout << "deltaEnergy/dt  " << mainSys.deltaTotEnergy/pars.dt << std::endl;
////            std::cout << "L2NormResidual  " << mainSys.L2NormResidual << std::endl;
//
//            std::cout << "\n" << std::endl;
//
//            if ( (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) > 1E10  || (mainSys.forceX.dot(mainSys.forceX)+ mainSys.forceY.dot(mainSys.forceY)) < 1E-10 || mainSys.maxR>50000.0){
//                std::cout << "Foce condition met !" << std::endl;
//                break;
//            }
//            if ( isnan(mainSys.areaRatio.sum()) || isnan(mainSys.forceX.sum()) ||  isnan(mainSys.forceY.minCoeff()) ||  isnan(mainSys.maxR) ){
//                std::cout << "System blew up !" << std::endl;
//                break;
//            }
//        }
    else if (pars.runMode=="continuousShear"){  //  shearing ***********************************************************************************************
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
        if (pars.boundaryType == "walls"){
            mainSys.compute_forces_harmonicWalls(baseData, pars, timeStep, 1, 0, pars.calculateHessian);
        }else if (pars.boundaryType == "periodic"){
            mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, pars.calculateHessian);
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
            mainSys.dump_global_data(pars, timeStep, "append", "running");


            if (timeStep % pars.dumpEvery == 0) {
                mainSys.dump_per_node(baseData, pars, timeStep);
                mainSys.dump_per_ele(baseData, pars,timeStep);
                if (pars.dumpPeriodicImagesXY){
                    mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
                }
                if (pars.identifyAndDumbFacets) {
                    mainSys.dump_facets(baseData, pars, timeStep);
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
