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
#include "integrators.hpp"

void fire2_solver(const BaseSysData& baseData, const Parameters& pars, long& timeStep , Configuration& mainSys){
    
    // simple parameters validation
    if (pars.FIRE_dtmax < pars.FIRE_dtmin) {
        std::cout << " FIRE_dtmax has to be larger than FIRE_dtmin" << std::endl;
        exit(1);
    }
    if (pars.FIRE_finc < 1.0) {
        std::cout << " FIRE_finc has to be greater than 1.0" << std::endl;
        exit(1);
    }
    
    if (pars.FIRE_fdec > 1.0) {
        std::cout << " FIRE_fdec has to be less than 1.0" << std::endl;
        exit(1);
    }
    
    // initialize some varialbles
    double FIRE_alpha = pars.FIRE_alpha_start;
    double FIRE_N_positive=0;
    double FIRE_N_negative=0;
    double FIRE_dt = pars.FIRE_dt_start;
    double FdotF;
    double VdotV;
    double scale1=0;
    double scale2=0;
    
    // initialize the system
    if (pars.boundaryType == "walls"){
        mainSys.compute_forces_harmonic_walls(baseData, pars, timeStep, 1, 0, pars.calculateHessian);
    }else if (pars.boundaryType == "periodic"){
        mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, pars.calculateHessian);
    }
    
    // initialize velocities
    mainSys.velocityX.fill(0);
    mainSys.velocityY.fill(0);
   
    // Leap Frog integration initialization
    if (pars.integrator == 1){
        mainSys.velocityX -= 0.5*FIRE_dt * mainSys.forceX;
        mainSys.velocityY -= 0.5*FIRE_dt * mainSys.forceY;
    }


    for (long i=1; i<= pars.FIRE_Nmax; i++) {
       

        if (pars.boundaryType == "walls"){
            mainSys.compute_forces_harmonic_walls(baseData, pars, timeStep, 1, 0, pars.calculateHessian);
        }else if (pars.boundaryType == "periodic"){
            mainSys.compute_forces_pbc(baseData, pars, timeStep, 1, 1, pars.calculateHessian);
        }
       
        
        mainSys.update_post_processing_data(baseData, pars);
        
        std::cout << "timeStep  " << timeStep << std::endl;
        std::cout << "e0  " << mainSys.e0 << std::endl;
        std::cout << "e1  " << mainSys.e1 << std::endl;
        std::cout << "phi  " << mainSys.phi << std::endl;
        std::cout << "energy  " <<  std::setprecision(12) << mainSys.totalEnergy <<std::endl;
        std::cout << "maxForce  " << mainSys.maxR << std::endl;
        std::cout << "avgForce  " << mainSys.avgR << std::endl;
        
        if (mainSys.maxR <= pars.RTolerance){
            std::cout << " Done !" << std::endl;
            break;
        }else if ( isnan(mainSys.areaRatio.sum()) || isnan(mainSys.forceX.sum()) ||  isnan(mainSys.forceY.minCoeff())){
            std::cout << "System blew up !" << std::endl;
            exit(1);
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
        
        //MD integration
        if (pars.integrator==0) semi_implicit_Euler(mainSys, FIRE_dt, 1, power, scale1, scale2);
        else if (pars.integrator==1) leapfrog(mainSys, FIRE_dt, 1, power, scale1, scale2);
        else if (pars.integrator==2) explicit_Euler(mainSys, FIRE_dt, 1, power, scale1, scale2);
        else if (pars.integrator==3) velocity_Verlet(mainSys, FIRE_dt, 1, power, scale1, scale2);
        
        
        std::cout << "power   " << power <<std::endl;
        std::cout << "FIRE_dt   " << FIRE_dt <<std::endl;
        std::cout << "FIRE_alpha   " << FIRE_alpha  <<std::endl;
        std::cout << "FIRE_Np+   " << FIRE_N_positive  <<std::endl;
        std::cout << "FIRE_Np-   " << FIRE_N_negative  <<std::endl;
        std::cout << "\n" << std::endl;
        
        timeStep++;
        
    }
        
}
