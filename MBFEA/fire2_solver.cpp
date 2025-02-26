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

void fire2_solver(const BaseSysData& baseData, const Parameters& pars, long& timeStep ,  Configuration& mainSys, long step, bool surfaceInteractions){
    
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
    long FIRE_N_positive=0;
    long FIRE_N_negative=0;
    double FIRE_dt = pars.FIRE_dt_start;
    double FdotF;
    double VdotV;
    double scale1=0;
    double scale2=0;
    double vmax=0;
    
    // initialize the system
    

    if (pars.boundaryType == "walls"){
        mainSys.compute_forces_walls(baseData, pars, timeStep, surfaceInteractions, 0, pars.calculateHessian);
    }else if (pars.boundaryType == "periodic"){
        mainSys.compute_forces_pbc(baseData, pars, timeStep, surfaceInteractions, 1, pars.calculateHessian);
    }
   
    // initialize velocities
    mainSys.velocityX.fill(0);
    mainSys.velocityY.fill(0);
   
    // Leap Frog integration initialization
    if (pars.integrator == 1){
        mainSys.velocityX -= 0.5*FIRE_dt * mainSys.forceX;
        mainSys.velocityY -= 0.5*FIRE_dt * mainSys.forceY;
        
        if (pars.runMode=="surfaceShear"){
               for (int nodeID=0; nodeID < baseData.numOriginalSurfaceNodes; nodeID++){
                   mainSys.velocityX[baseData.flatSurfaceNodes[nodeID]]=0;
                   mainSys.velocityY[baseData.flatSurfaceNodes[nodeID]]=0;
               }
        }
    }


    for (long i=1; i<= pars.FIRE_Nmax; i++) {
        

        if (pars.boundaryType == "walls"){
            mainSys.compute_forces_walls(baseData, pars, timeStep, surfaceInteractions, 0, pars.calculateHessian);
        }else if (pars.boundaryType == "periodic"){
            mainSys.compute_forces_pbc(baseData, pars, timeStep, surfaceInteractions, 1, pars.calculateHessian);
        }
       
        if (pars.runMode=="surfaceShear"){
               for (int nodeID=0; nodeID < baseData.numSurfaceNodes; nodeID++){
                   mainSys.forceX[baseData.flatSurfaceNodes[nodeID]]=0;
                   mainSys.forceY[baseData.flatSurfaceNodes[nodeID]]=0;
               }
        }
        
        mainSys.update_post_processing_data(baseData, pars);
        
        
        
        if (mainSys.maxR <= pars.FIRE_RTolerance){
            std::cout << "step " << step << " is done !" << std::endl;
           
            break;
        }else if ( isnan(mainSys.areaRatio.sum()) || isnan(mainSys.forceX.sum()) ||  isnan(mainSys.forceY.minCoeff())){
            std::cout << "System blew up !" << std::endl;
            std::cout << "mainSys.areaRatio.sum() = " << mainSys.areaRatio.sum()<<std::endl;
            std::cout << "SmainSys.forceX.sum() = " << mainSys.forceX.sum() << std::endl;
            std::cout << "mainSys.forceY.minCoeff()  = " << mainSys.forceY.minCoeff() << std::endl;
            exit(1);
        }
        
        double power = mainSys.forceX.dot(mainSys.velocityX) + mainSys.forceY.dot(mainSys.velocityY);
        
        if(timeStep%pars.writeToConsoleEvery==0){
            
            std::cout << "step  " << step << "\n" << std::endl;
            std::cout << "timeStep  " << timeStep << std::endl;
            std::cout << "energy  " <<  std::setprecision(12) << mainSys.totalEnergy <<std::endl;
            std::cout << "power   " << power <<std::endl;
            std::cout << "e0  " << mainSys.e0_W << std::endl;
            std::cout << "e1  " << mainSys.e1_W << std::endl;
            std::cout << "phi  " << mainSys.phi << std::endl;
            std::cout << "pressure  " << mainSys.P_BB << std::endl;
            std::cout << "maxForce  " << mainSys.maxR << std::endl;
            std::cout << "avgForce  " << mainSys.avgR << std::endl;
            std::cout << "L2R  " << mainSys.L2NormResidual << std::endl;
            std::cout << "interactions  " << mainSys.nodeIinteractions + mainSys.segmentIinteractions << std::endl;
            std::cout << "maximum peneteration  " << mainSys.maxInterference << std::endl;
            std::cout << "areaRatio.sum() = " << mainSys.areaRatio.sum()<<std::endl;
            std::cout << "forceX.sum() = " << mainSys.forceX.sum() << std::endl;
            std::cout << "forceY.minCoeff()  = " << mainSys.forceY.minCoeff() << std::endl;
            std::cout << "FIRE_dt   " << FIRE_dt <<std::endl;
            std::cout << "FIRE_alpha   " << FIRE_alpha  <<std::endl;
            std::cout << "FIRE_Np+   " << FIRE_N_positive  <<std::endl;
            std::cout << "FIRE_Np-   " << FIRE_N_negative  <<std::endl;
            std::cout << "\n" << std::endl;
        }
        
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
        vmax= sqrt((mainSys.velocityX.array()*mainSys.velocityX.array()+mainSys.velocityY.array()*mainSys.velocityY.array()).maxCoeff());
        
        while (vmax*FIRE_dt>pars.verletCellCutoff*1.5){
            if ( FIRE_dt*pars.FIRE_fdec >= pars.FIRE_dtmin){
                FIRE_dt*=pars.FIRE_fdec;
                std::cout << "FIRE_dt   " << FIRE_dt <<std::endl;
            }
        }
        //MD integration
        if (pars.integrator==0) semi_implicit_Euler(mainSys, FIRE_dt, 1, power, scale1, scale2);
        else if (pars.integrator==1) leapfrog(mainSys, FIRE_dt, 1, power, scale1, scale2);
        else if (pars.integrator==2) explicit_Euler(mainSys, FIRE_dt, 1, power, scale1, scale2);
        else if (pars.integrator==3) velocity_Verlet(mainSys, FIRE_dt, 1, power, scale1, scale2);
        
        
        timeStep++;
        
//        if(timeStep==139){
//             mainSys.dump_smoothcurves(baseData,pars,timeStep);
//             mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
//            mainSys.dump_per_node(baseData, pars, timeStep);
//            mainSys.dump_per_ele(baseData, pars,timeStep);
//
//        }
    }
        
}

