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

void gd_solver(const BaseSysData& baseData, const Parameters& pars, long& timeStep, Configuration& mainSys){

            // calculate forces and energy of the current configuration of the system
            if (pars.boundaryType == "walls"){
                mainSys.compute_forces_harmonic_walls(baseData, pars, timeStep, 1, 0, pars.calculateHessian);
            } else if (pars.boundaryType == "periodic"){
                mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, pars.calculateHessian);
            }
            // Postporcesseing calculations
            mainSys.update_post_processing_data(baseData, pars);
            //dump
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
            std::cout << "force tolerance  " << pars.maxForceTol << std::endl;
            std::cout << "maxForce  " << mainSys.maxR << std::endl;
            std::cout << "meanForce  " << mainSys.avgR << std::endl;
            std::cout << "phi  " << mainSys.phi << std::endl;
            std::cout << "e0  " << mainSys.e0 << std::endl;
            std::cout << "e1  " << mainSys.e1 << std::endl;
            std::cout << "\n" << std::endl;
    
            if ( mainSys.maxR<=pars.maxForceTol){
                std::cout << "Foce condition met !" << std::endl;
                
                mainSys.dump_per_node(baseData, pars, timeStep);
                mainSys.dump_per_ele(baseData, pars,timeStep);
                if (pars.dumpPeriodicImagesXY){
                    mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
                }
                if (pars.identifyAndDumbFacets) {
                    mainSys.dump_facets(baseData, pars, timeStep);
                }
                
            } else if ( isnan(mainSys.areaRatio.sum()) || isnan(mainSys.forceX.sum()) ||  isnan(mainSys.forceY.minCoeff())){
                std::cout << "System blew up !" << std::endl;
                exit(1);
            }
            
            // Take explicit Euler step
            explicit_Euler(mainSys, pars.dt, 0, 0, 0, 0);
        
            timeStep++;
}


