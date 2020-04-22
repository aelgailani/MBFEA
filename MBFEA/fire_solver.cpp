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

void fire_solver(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys){
    double FIRE_alpha =  pars.FIRE_alpha_start;
    double FIRE_N = 0;
    double FIRE_dt = pars.FIRE_dt_start;
    double FIRE_prevDt = pars.FIRE_dt_start;
    
    if (pars.runMode=="compress"){
        assert(pars.startingStrainStep <= pars.numStrainSteps);
        
        double target_e0 = - log(sqrt(pars.Ap/(pars.targetPhi*baseData.lxRef*baseData.lyRef))); // my conviension is positive e0 for compression

        for (long strainStep = pars.startingStrainStep ; strainStep<= pars.numStrainSteps ; strainStep++) {
                        
            mainSys.compress(baseData, pars, target_e0 * float(strainStep)/float(pars.numStrainSteps));

            while (1)
            {
                

                if (pars.boundaryType == "walls"){
                    mainSys.compute_forces_harmonic_walls(baseData, pars, timeStep, 1, 0, pars.calculateHessian);
                }else if (pars.boundaryType == "periodic"){
                    mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, pars.calculateHessian);
                }
                mainSys.update_post_processing_data(baseData, pars);
                
                std::cout << "timeStep  " << timeStep << std::endl;
                std::cout << "strainStep  " << strainStep<<"   out of  " << pars.numStrainSteps <<std::endl;
                std::cout << "e0  " << mainSys.e0 << ", target is  " << target_e0 <<std::endl;
                std::cout << "phi  " << mainSys.phi << ", target is  " << pars.targetPhi <<std::endl;
                std::cout << "energy  " <<  std::setprecision(12) << mainSys.totalEnergy <<std::endl;
                std::cout << "maxForce  " << mainSys.maxR << std::endl;
                std::cout << "avgForce  " << mainSys.avgR << std::endl;
                
                if ( mainSys.maxR > 1E10  || mainSys.maxR <= pars.RTolerance){
                    std::cout << "Foce condition met !" << std::endl;
                    break;
                }
                
                double power = mainSys.forceX.dot(mainSys.velocityX)+ mainSys.forceY.dot(mainSys.velocityY);

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
                
                
                mainSys.curPosX += mainSys.velocityX*FIRE_dt;
                mainSys.curPosY += mainSys.velocityY*FIRE_dt;
                mainSys.velocityX += mainSys.forceX*FIRE_dt;
                mainSys.velocityY += mainSys.forceY*FIRE_dt;
                
                std::cout << "power   " << power <<std::endl;
                std::cout << "FIRE_dt   " << FIRE_dt <<std::endl;
                std::cout << "FIRE_alpha   " << FIRE_alpha  <<std::endl;
                
                std::cout << "\n" << std::endl;
                timeStep++;
                
                

            }
            FIRE_alpha = pars.FIRE_alpha_start;
            FIRE_N = timeStep;
            mainSys.dump_global_data(pars, timeStep, "append", "final");
            mainSys.dump_per_node(baseData, pars, strainStep);
            mainSys.dump_per_ele(baseData, pars,strainStep);
            if (pars.dumpPeriodicImagesXY){
                mainSys.dump_per_node_periodic_images_on(baseData, pars, strainStep);
            }
           
        }
        
    }else if (pars.runMode=="stepShear"){
        
        long stage=0;
        mainSys.dump_global_data(pars, timeStep, "write", "final"); //open file and write cols names
        while(stage<2){
            
            mainSys.shear(baseData, pars, pars.targetShear);
            
            double FIRE_alpha =  pars.FIRE_alpha_start;
            double FIRE_N = 0;
            double FIRE_dt = pars.FIRE_dt_start;
            double FIRE_prevDt = pars.FIRE_dt_start;
            
            while (1)
            {
                
                if (pars.boundaryType == "walls"){
                    mainSys.compute_forces_harmonic_walls(baseData, pars, timeStep, 1, 0, pars.calculateHessian);
                }else if (pars.boundaryType == "periodic"){
                    mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, pars.calculateHessian);
                }
                mainSys.update_post_processing_data(baseData, pars);
                
                std::cout << "timeStep  " << timeStep << std::endl;
                std::cout << "stage  " << stage << std::endl;
                std::cout << "e1  " << mainSys.e1 << std::endl;
                std::cout << "phi  " << mainSys.phi <<std::endl;
                std::cout << "energy  " <<  std::setprecision(12) << mainSys.totalEnergy <<std::endl;
                std::cout << "maxForce  " << mainSys.maxR << std::endl;
                std::cout << "avgForce  " << mainSys.avgR << std::endl;
                
                if ( mainSys.maxR > 1E10  || mainSys.maxR <= pars.RTolerance){
                   std::cout << "Foce condition met !" << std::endl;
                   break;
                }
                
                double power = mainSys.forceX.dot(mainSys.velocityX)+ mainSys.forceY.dot(mainSys.velocityY);

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
                
                
               
                mainSys.curPosX += mainSys.velocityX*FIRE_dt;
                mainSys.curPosY += mainSys.velocityY*FIRE_dt;
                mainSys.velocityX += mainSys.forceX*FIRE_dt;
                mainSys.velocityY += mainSys.forceY*FIRE_dt;
                
                std::cout << "power   " << power <<std::endl;
                std::cout << "FIRE_dt   " << FIRE_dt <<std::endl;
                std::cout << "FIRE_alpha   " << FIRE_alpha  <<std::endl;
                std::cout << "\n" << std::endl;
                
                timeStep++;
            }
            FIRE_alpha = pars.FIRE_alpha_start;
            FIRE_N = timeStep;
            mainSys.dump_global_data(pars, timeStep, "append", "final");
            mainSys.dump_per_node(baseData, pars, stage);
            mainSys.dump_per_ele(baseData, pars,stage);
            if (pars.dumpPeriodicImagesXY){
                mainSys.dump_per_node_periodic_images_on(baseData, pars, stage);
            }
            
            stage+=1;
            
        }
    
    }
        
}
