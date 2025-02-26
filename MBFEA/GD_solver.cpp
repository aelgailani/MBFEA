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

void gd_solver(const BaseSysData& baseData, const Parameters& pars, long& timeStep, long& step, std::string name, Configuration& mainSys, bool dumpStateData, bool surfaceInteractions){
            
            // calculate forces and energy of the current configuration of the system
            if (pars.boundaryType == "walls"){
                mainSys.compute_forces_walls(baseData, pars, timeStep, surfaceInteractions, 0, pars.calculateHessian);
            } else if (pars.boundaryType == "periodic"){
                mainSys.compute_forces_pbc(baseData, pars, timeStep, surfaceInteractions, 1, pars.calculateHessian);
            }
    
            // zero out surface nodes forces if you want them fixed in their homo position
            if (pars.runMode=="surfaceShear"){
                for (int nodeID=0; nodeID < baseData.numOriginalSurfaceNodes; nodeID++){
                    mainSys.forceX[baseData.flatSurfaceNodes[nodeID]]=0;
                    mainSys.forceY[baseData.flatSurfaceNodes[nodeID]]=0;
                }
            }
            // Postporcesseing calculations
            mainSys.update_post_processing_data(baseData, pars);
            //dump
            mainSys.dump_global_data(pars, timeStep, name, "append", "running");
            
            if (dumpStateData==true && timeStep % pars.dumpEvery == 0) {
                mainSys.dump_per_node(baseData, pars, step);
                mainSys.dump_per_ele(baseData, pars,step);
                if (pars.dumpPeriodicImagesXY){
                    mainSys.dump_per_node_periodic_images_on(baseData, pars, step);
                }
                if (pars.identifyAndDumbFacets) {
                    if (pars.contactMethod=="nts"){
                        mainSys.dump_facets(baseData, pars, step);
                    }else{
                        mainSys.dump_facets_ntn(baseData, pars, step);
                    }
                    
                }
                if(pars.dumpSmoothenCurves){
                    mainSys.dump_smoothcurves(baseData,pars,step);
                }
                step++;
            }
//            if(timeStep%pars.writeToConsoleEvery==0){
//                std::cout << "force tolerance  " << pars.maxForceTol << std::endl;
//                std::cout << "maxForce  " << mainSys.maxR << std::endl;
//                std::cout << "meanForce  " << mainSys.avgR << std::endl;
//                std::cout << "L2NormR  " << mainSys.L2NormResidual << std::endl;
//                std::cout << "\n" << std::endl;
//            }
            if ( mainSys.maxR<=pars.maxForceTol){
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


