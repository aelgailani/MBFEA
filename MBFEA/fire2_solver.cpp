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

void fire2_solver(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys){
    double FIRE_alpha =  pars.FIRE_alpha_start;
    double FIRE_N_positive = 0; //steps since positve power
    double FIRE_N_negative = 0;//steps since negative power
    double FIRE_dt = pars.FIRE_dt_start;
    mainSys.velocityX.fill(0);
    mainSys.velocityY.fill(0);
    double FdotF;
    double VdotV;
    double scale1;
    double scale2;
    if (pars.runMode=="compress"){
        assert(pars.startingStrainStep <= pars.numStrainSteps);
        double refPhi = pars.Ap/(baseData.lxRef*baseData.lyRef);
        double target_e0 = - log(sqrt(pars.Ap/(pars.targetPhi*baseData.lxRef*baseData.lyRef))); // my conviension is positive e0 for compression

        for (long strainStep = pars.startingStrainStep ; strainStep<= pars.numStrainSteps ; strainStep++) {
            
            double stepPhi = refPhi + (pars.targetPhi - refPhi) * float(strainStep)/float(pars.numStrainSteps); // the required phi of this step as a fraction of the required target phi
            double relativeStrain = - log(sqrt(pars.Ap/(stepPhi*baseData.lxRef*baseData.lyRef))) - mainSys.e0;  // strain between current and new configuration
            mainSys.compress(baseData, pars, relativeStrain);

            for (long i=1; i<= pars.FIRE_Nmax; i++) {
               

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
                
                if (mainSys.maxR <= pars.RTolerance){
                    std::cout << " Done !" << std::endl;
                    break;
                }
                
                double power = mainSys.forceX.dot(mainSys.velocityX) + mainSys.forceY.dot(mainSys.velocityY);

                if (power >0) {
                    FIRE_N_positive +=1;
                    FIRE_N_negative =0;
                    if (FIRE_N_positive > pars.FIRE_N_positive_min){
                        FIRE_dt = fmin(FIRE_dt * pars.FIRE_finc, pars.FIRE_dtmax);
                        FIRE_alpha *= pars.FIRE_falpha;
                    }
                    scale1=(1-FIRE_alpha);
                    FdotF = mainSys.forceX.dot(mainSys.forceX)+mainSys.forceY.dot(mainSys.forceY);
                    VdotV = mainSys.velocityX.dot(mainSys.velocityX)+mainSys.velocityY.dot(mainSys.velocityY);
                    if (FdotF <= 1e-20) scale2 = 0.0;
                    else scale2 = FIRE_alpha * sqrt(VdotV/FdotF);
                    
                } else {
                    
                    FIRE_N_positive =0;
                    FIRE_N_negative +=1;
                    
                    if (FIRE_N_negative > pars.FIRE_N_negative_max){
                        std::cout << " Failed to converge !" << std::endl;
                        break;
                    }
                    if (!(pars.FIRE_intialdelay && i < pars.FIRE_N_positive_min)){
                        if ( FIRE_dt*pars.FIRE_fdec >= pars.FIRE_dtmin){
                            FIRE_dt*=pars.FIRE_fdec;
                        }
                        FIRE_alpha = pars.FIRE_alpha_start;
                    }
                    
                    mainSys.curPosX -= 0.5* mainSys.velocityX*FIRE_dt;
                    mainSys.curPosY -= 0.5* mainSys.velocityY*FIRE_dt;
                    mainSys.velocityX.fill(0);
                    mainSys.velocityY.fill(0);
                    
                }
                
                // Semi-implicit Euler OR Leap Frog integration
                mainSys.velocityX += mainSys.forceX*FIRE_dt;
                mainSys.velocityY += mainSys.forceY*FIRE_dt;
                
                if (power>0.0){
                    mainSys.velocityX = scale1 * mainSys.velocityX+ scale2 * mainSys.forceX;
                    mainSys.velocityY = scale1 * mainSys.velocityY+ scale2 * mainSys.forceY;
                }
                mainSys.curPosX += mainSys.velocityX*FIRE_dt;
                mainSys.curPosY += mainSys.velocityY*FIRE_dt;
                
                
                std::cout << "power   " << power <<std::endl;
                std::cout << "FIRE_dt   " << FIRE_dt <<std::endl;
                std::cout << "FIRE_alpha   " << FIRE_alpha  <<std::endl;
                std::cout << "FIRE_Np+   " << FIRE_N_positive  <<std::endl;
                std::cout << "FIRE_Np-   " << FIRE_N_negative  <<std::endl;
                std::cout << "\n" << std::endl;
                timeStep++;
                
                

            }
            
            if (mainSys.maxR > pars.RTolerance){
                std::cout << " Failed to coverge !" << std::endl;

            }else{
                mainSys.dump_global_data(pars, timeStep, "append", "final");
                mainSys.dump_per_node(baseData, pars, strainStep);
                mainSys.dump_per_ele(baseData, pars,strainStep);
                if (pars.dumpPeriodicImagesXY){
                    mainSys.dump_per_node_periodic_images_on(baseData, pars, strainStep);
                }
            }
            
           
            FIRE_alpha = pars.FIRE_alpha_start;
            FIRE_N_positive=0;
            FIRE_N_negative=0;
            FIRE_dt = pars.FIRE_dt_start;
            mainSys.velocityX.fill(0);
            mainSys.velocityY.fill(0);

        }
    }
    else if (pars.runMode=="stepShear"){

        long stage=0;
        mainSys.dump_global_data(pars, timeStep, "write", "final"); //open file and write cols names
        while(stage<2){

            mainSys.shear(baseData, pars, pars.targetShear);


            for (long i=1; i<= pars.FIRE_Nmax; i++) {
               

                if (pars.boundaryType == "walls"){
                    mainSys.compute_forces_harmonic_walls(baseData, pars, timeStep, 1, 0, pars.calculateHessian);
                }else if (pars.boundaryType == "periodic"){
                    mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, pars.calculateHessian);
                }
               
                
                mainSys.update_post_processing_data(baseData, pars);
                
                std::cout << "timeStep  " << timeStep << std::endl;
                std::cout << "stage  " << stage <<std::endl;
                std::cout << "e1  " << mainSys.e1 << ", target is  " << pars.targetShear*2 <<std::endl;
                std::cout << "phi  " << mainSys.phi <<std::endl;
                std::cout << "energy  " <<  std::setprecision(12) << mainSys.totalEnergy <<std::endl;
                std::cout << "maxForce  " << mainSys.maxR << std::endl;
                std::cout << "avgForce  " << mainSys.avgR << std::endl;
                
                if (mainSys.maxR <= pars.RTolerance){
                    std::cout << " Done !" << std::endl;
                    break;
                }
                
                double power = mainSys.forceX.dot(mainSys.velocityX) + mainSys.forceY.dot(mainSys.velocityY);

                if (power >0) {
                    FIRE_N_positive +=1;
                    FIRE_N_negative =0;
                    if (FIRE_N_positive > pars.FIRE_N_positive_min){
                        FIRE_dt = fmin(FIRE_dt * pars.FIRE_finc, pars.FIRE_dtmax);
                        FIRE_alpha *= pars.FIRE_falpha;
                    }
                    scale1=(1-FIRE_alpha);
                    FdotF = mainSys.forceX.dot(mainSys.forceX)+mainSys.forceY.dot(mainSys.forceY);
                    VdotV = mainSys.velocityX.dot(mainSys.velocityX)+mainSys.velocityY.dot(mainSys.velocityY);
                    if (FdotF <= 1e-20) scale2 = 0.0;
                    else scale2 = FIRE_alpha * sqrt(VdotV/FdotF);
                    
                } else {
                    
                    FIRE_N_positive =0;
                    FIRE_N_negative +=1;
                    
                    if (FIRE_N_negative > pars.FIRE_N_negative_max){
                        std::cout << " Failed to converge !" << std::endl;
                        break;
                    }
                    if (!(pars.FIRE_intialdelay && i < pars.FIRE_N_positive_min)){
                        if ( FIRE_dt*pars.FIRE_fdec >= pars.FIRE_dtmin){
                            FIRE_dt*=pars.FIRE_fdec;
                        }
                        FIRE_alpha = pars.FIRE_alpha_start;
                    }
                    
                    mainSys.curPosX -= 0.5* mainSys.velocityX*FIRE_dt;
                    mainSys.curPosY -= 0.5* mainSys.velocityY*FIRE_dt;
                    mainSys.velocityX.fill(0);
                    mainSys.velocityY.fill(0);
                    
                }
                
                // Semi-implicit Euler OR Leap Frog integration
                mainSys.velocityX += mainSys.forceX*FIRE_dt;
                mainSys.velocityY += mainSys.forceY*FIRE_dt;
                
                if (power>0.0){
                    mainSys.velocityX = scale1 * mainSys.velocityX+ scale2 * mainSys.forceX;
                    mainSys.velocityY = scale1 * mainSys.velocityY+ scale2 * mainSys.forceY;
                }
                mainSys.curPosX += mainSys.velocityX*FIRE_dt;
                mainSys.curPosY += mainSys.velocityY*FIRE_dt;
                
                
                std::cout << "power   " << power <<std::endl;
                std::cout << "FIRE_dt   " << FIRE_dt <<std::endl;
                std::cout << "FIRE_alpha   " << FIRE_alpha  <<std::endl;
                std::cout << "FIRE_Np+   " << FIRE_N_positive  <<std::endl;
                std::cout << "FIRE_Np-   " << FIRE_N_negative  <<std::endl;
                std::cout << "\n" << std::endl;
                timeStep++;
                
                

            }
            
           if (mainSys.maxR > pars.RTolerance){
                 std::cout << " Failed to coverge !" << std::endl;

             }else{
                 mainSys.dump_global_data(pars, timeStep, "append", "final");
                 mainSys.dump_per_node(baseData, pars, stage);
                 mainSys.dump_per_ele(baseData, pars,stage);
                 if (pars.dumpPeriodicImagesXY){
                     mainSys.dump_per_node_periodic_images_on(baseData, pars, stage);
                 }
             }
             
            
             FIRE_alpha = pars.FIRE_alpha_start;
             FIRE_N_positive=0;
             FIRE_N_negative=0;
             FIRE_dt = pars.FIRE_dt_start;
             mainSys.velocityX.fill(0);
             mainSys.velocityY.fill(0);

            stage+=1;

        }

    }
        
}
