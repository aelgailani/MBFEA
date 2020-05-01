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
    double FIRE_Np = 0; //stepsSincePositvePower
    double FIRE_dt = pars.FIRE_dt_start;
    double FIRE_prevDt = pars.FIRE_dt_start;
    double force_magnitude;
    double velocity_magnitude;
    double scale1;
    double scale2;
    if (pars.runMode=="compress"){
        assert(pars.FIRE_startingStrainStep <= pars.FIRE_numStrainSteps);
        double refPhi = pars.Ap/(baseData.lxRef*baseData.lyRef);
        double target_e0 = - log(sqrt(pars.Ap/(pars.targetPhi*baseData.lxRef*baseData.lyRef))); // my conviension is positive e0 for compression

        for (long strainStep = pars.FIRE_startingStrainStep ; strainStep<= pars.FIRE_numStrainSteps ; strainStep++) {
            
            double stepPhi = refPhi + (pars.targetPhi - refPhi) * float(strainStep)/float(pars.FIRE_numStrainSteps); // the required phi of this step as a fraction of the required target phi
            double relativeStrain = - log(sqrt(pars.Ap/(stepPhi*baseData.lxRef*baseData.lyRef))) - mainSys.e0;  // strain between current and new configuration
            mainSys.affine_compression(baseData, pars, relativeStrain);

            while (1)
            {
                

                if (pars.boundaryType == "walls"){
                    mainSys.compute_forces_harmonic_walls(baseData, pars, timeStep, 1, 0, pars.calculateHessian);
                }else if (pars.boundaryType == "periodic"){
                    mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, pars.calculateHessian);
                }
               
                
                mainSys.update_post_processing_data(baseData, pars);
                
                std::cout << "timeStep  " << timeStep << std::endl;
                std::cout << "strainStep  " << strainStep<<"   out of  " << pars.FIRE_numStrainSteps <<std::endl;
                std::cout << "e0  " << mainSys.e0 << ", target is  " << target_e0 <<std::endl;
                std::cout << "phi  " << mainSys.phi << ", target is  " << pars.targetPhi <<std::endl;
                std::cout << "energy  " <<  std::setprecision(12) << mainSys.totalEnergy <<std::endl;
                std::cout << "maxForce  " << mainSys.maxR << std::endl;
                std::cout << "avgForce  " << mainSys.avgR << std::endl;
                
                if ( mainSys.maxR > 1E10  || mainSys.maxR <= pars.FIRE_RTolerance){
                    std::cout << "Foce condition met !" << std::endl;
                    break;
                }
                
                double power = mainSys.forceX.dot(mainSys.velocityX)+ mainSys.forceY.dot(mainSys.velocityY);

                if (power <=0) {
                    FIRE_dt *= pars.FIRE_fdec;
                    mainSys.velocityX.fill(0);
                    mainSys.velocityY.fill(0);
                    FIRE_alpha = pars.FIRE_alpha_start;
                    FIRE_Np=0;
                } else {
                    FIRE_Np +=1;
                    if (FIRE_Np > pars.FIRE_N_positive_min){
                        FIRE_prevDt = FIRE_dt;
                        FIRE_dt = fmin(FIRE_dt * pars.FIRE_finc, pars.FIRE_dtmax);
                        FIRE_alpha *= pars.FIRE_falpha;
                    }
                    scale1=(1-FIRE_alpha);
                    force_magnitude = (mainSys.forceX.array().pow(2)+mainSys.forceY.array().pow(2)).sum();
                    velocity_magnitude = (mainSys.velocityX.array().pow(2)+mainSys.velocityY.array().pow(2)).sum();
                    if (force_magnitude <= 1e-20) scale2 = 0.0;
                    else scale2 = FIRE_alpha * sqrt(velocity_magnitude/force_magnitude);
                    mainSys.velocityX = scale1 * mainSys.velocityX+ scale2 * mainSys.forceX;
                    mainSys.velocityY = scale1 * mainSys.velocityY+ scale2 * mainSys.forceY;
                }
                
                mainSys.velocityX += mainSys.forceX*FIRE_dt;
                mainSys.velocityY += mainSys.forceY*FIRE_dt;
                mainSys.curPosX += mainSys.velocityX*FIRE_dt;
                mainSys.curPosY += mainSys.velocityY*FIRE_dt;
                
                
                
                std::cout << "power   " << power <<std::endl;
                std::cout << "FIRE_dt   " << FIRE_dt <<std::endl;
                std::cout << "FIRE_alpha   " << FIRE_alpha  <<std::endl;
                std::cout << "FIRE_Np   " << FIRE_Np  <<std::endl;
                
                std::cout << "\n" << std::endl;
                timeStep++;
                
                

            }
            FIRE_alpha = pars.FIRE_alpha_start;
            FIRE_Np=0;
            FIRE_dt = pars.FIRE_dt_start;
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
            
            mainSys.affine_axial_shearing(baseData, pars, pars.targetShear);
            
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
                
                if ( mainSys.maxR > 1E10  || mainSys.maxR <= pars.FIRE_RTolerance){
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
                    if ((timeStep - FIRE_N) > pars.FIRE_N_positive_min){
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
