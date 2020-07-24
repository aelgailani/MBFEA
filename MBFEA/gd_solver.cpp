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
#include "integrators.hpp"

void gd_solver(const BaseSysData& baseData, const Parameters& pars, long& timeStep, std::string name, Configuration& mainSys, bool dumpStateData){

            // calculate forces and energy of the current configuration of the system
            if (pars.boundaryType == "walls"){
                mainSys.compute_forces_walls(baseData, pars, timeStep, 1, 0, pars.calculateHessian);
            } else if (pars.boundaryType == "periodic"){
                mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, pars.calculateHessian);
            }
            // Postporcesseing calculations
            mainSys.update_post_processing_data(baseData, pars);
            //dump
            mainSys.dump_global_data(pars, timeStep, name, "append", "running");
            
            if (dumpStateData==true && timeStep % pars.dumpEvery == 0) {
                mainSys.dump_per_node(baseData, pars, timeStep);
                mainSys.dump_per_ele(baseData, pars,timeStep);
                if (pars.dumpPeriodicImagesXY){
                    mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
                }
                if (pars.identifyAndDumbFacets) {
                    if (pars.contactMethod=="nts"){
                        mainSys.dump_facets(baseData, pars, timeStep);
                    }else{
                        mainSys.dump_facets_ntn(baseData, pars, timeStep);
                    }
                    
                }
                if(pars.dumpSmoothenCurves){
                    mainSys.dump_smoothcurves(baseData,pars,timeStep);
                }

            }
//            if(timeStep%pars.writeToConsoleEvery==0){
//                std::cout << "force tolerance  " << pars.maxForceTol << std::endl;
//                std::cout << "maxForce  " << mainSys.maxR << std::endl;
//                std::cout << "meanForce  " << mainSys.avgR << std::endl;
//                std::cout << "L2NormR  " << mainSys.L2NormResidual << std::endl;
//                std::cout << "\n" << std::endl;
//            }
            if ( mainSys.L2NormResidual<=pars.maxForceTol){
                std::cout << "Foce condition met !" << std::endl;
                
                
            } else if ( isnan(mainSys.areaRatio.sum()) || isnan(mainSys.forceX.sum()) ||  isnan(mainSys.forceY.minCoeff())){
    
                std::cout << "System blew up !" << std::endl;
                std::cout << "maxForce  " << mainSys.maxR << std::endl;
                std::cout << "meanForce  " << mainSys.avgR << std::endl;
                std::cout << "L2NormR  " << mainSys.L2NormResidual << std::endl;
                exit(1);
            }
            
            // Take explicit Euler step
        mainSys.curPosX += mainSys.forceX * pars.dt;
        mainSys.curPosY += mainSys.forceY *  pars.dt;
        
        timeStep++;
}


