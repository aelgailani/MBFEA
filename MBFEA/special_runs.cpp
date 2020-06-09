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
#include "standard_runs.hpp"

bool isInsideTriangle(float x, float xo, float x1, float x2, float y, float yo,float y1,float y2){
    float v1=((x-xo)*(yo-y1)+(y-yo)*(x1-xo));
    float v2=((x-x1)*(y1-y2)+(y-y1)*(x2-x1));
    float v3=((x-x2)*(y2-yo)+(y-y2)*(xo-x2));
    
    if (signbit(v1)==signbit(v2)==signbit(v3)) return true;
        else return false;
    
}

struct PowerlawRepulsion compute_powerlaw_replusion_by_segment(double x, double x1, double x2, double y, double y1,double y2, double sigma, double epsilon, double rcut){
    struct PowerlawRepulsion F;
    double Nx = (y2-y1);
    double Ny = - (x2-x1);
    double N_norm = sqrt(Nx*Nx+Ny*Ny);
    double nx = Nx/N_norm;
    double ny = Ny/N_norm;
    double r = (x-x1)*nx+(y-y1)*ny;
    if (r<rcut){
        double f = 0.5*epsilon/sigma*12*pow((sigma/r),13);
        F.fx = f * (nx);
        F.fy = f * (ny);
        F.energy = 0.5*epsilon*pow((sigma/r),12);
//        std::cout<< "r="<< r<< std::endl;
//        std::cout<< "Nx="<< Nx<< std::endl;
//        std::cout<< "Ny="<< Ny<< std::endl;
//        std::cout<< "fx="<< F.fx<< std::endl;
//        std::cout<< "fy="<< F.fy<< std::endl;
    }
    
    return F;
}

void shear_special_FIRE(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys){
    float shearStep=pars.deformationRate;
    long numSteps=100000;
    long step = 0;
    double refPhi = pars.Ap/(0.5*baseData.lxRef*baseData.lyRef);
    double target_e0 = - log(sqrt(refPhi/(pars.targetPhi)));
    double stepPhi = refPhi;
    // initial compression
    for (long strainStep = pars.FIRE_startingStrainStep ; strainStep<= pars.FIRE_numStrainSteps ; strainStep++) {
        double prevStepPhi = stepPhi;
        stepPhi = refPhi + (pars.targetPhi - refPhi) * float(strainStep)/float(pars.FIRE_numStrainSteps); // the required phi of this step as a fraction of the required target phi
        double relativeStrain = log(sqrt(prevStepPhi/stepPhi));  // strain between current and new configuration
        
        std::cout << timeStep << std::endl;
        std::cout << "phi  " << mainSys.phi * 2 << "  target is  " << pars.targetPhi << std::endl;
        std::cout << "e0  " << mainSys.e0 << "  target is  " << target_e0 << std::endl;
        
        mainSys.affine_compression_triWalls(baseData, pars, relativeStrain,1.75*pars.initialStretch, 1*pars.initialStretch);
        
        fire2_solver(baseData, pars, timeStep , mainSys, strainStep);
        
        if (mainSys.maxR > pars.FIRE_RTolerance){
             std::cout << " Failed to coverge compression !" << std::endl;

         }
//            else{
//             mainSys.dump_global_data(pars, timeStep, "append", "final");
//             mainSys.dump_per_node(baseData, pars, strainStep);
//             mainSys.dump_per_ele(baseData, pars,strainStep);
//             if (pars.dumpPeriodicImagesXY){
//                 mainSys.dump_per_node_periodic_images_on(baseData, pars, strainStep);
//             }
//         }
        if (mainSys.P2<pars.targetPressure) break;
    }
    mainSys.dump_global_data(pars, timeStep, "write", "final");
    while(step<=numSteps && mainSys.e1 <= pars.shearTo){
        
        
        
        
        std::cout << timeStep << std::endl;
        std::cout << "phi  " << mainSys.phi *2<< std::endl;
        std::cout << "e1  " << mainSys.e1 << "  target is  " << target_e0 << std::endl;
        
        mainSys.affine_axial_shearing_triWalls(baseData, pars, shearStep,1.75*pars.initialStretch, 1*pars.initialStretch);

        
        fire2_solver(baseData, pars, timeStep , mainSys, step);
        
        if (mainSys.maxR > pars.FIRE_RTolerance){
             std::cout << " Failed to coverge !" << std::endl;
             exit(1);
         }else{
             mainSys.dump_global_data(pars, timeStep, "append", "final");
             mainSys.dump_per_node(baseData, pars, step);
             mainSys.dump_per_ele(baseData, pars,step);
             if (pars.dumpPeriodicImagesXY){
                 mainSys.dump_per_node_periodic_images_on(baseData, pars, step);
             }
             
  
        }
        
        step++;
    }
    
}

void shear_special_GD(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys){

    double refPhi = pars.Ap/(0.5*baseData.lxRef*baseData.lyRef);
    double target_e0 = - log(sqrt(refPhi/(pars.targetPhi)));

    // initial compression
    gd_solver(baseData,pars,timeStep, mainSys, false);

    while(mainSys.phi <= pars.targetPhi){
        
        std::cout << timeStep << std::endl;
        std::cout << "phi   " << mainSys.phi << "  target is  " << pars.targetPhi << std::endl;
        std::cout << "e0   " << mainSys.e0 << "  target is  " << target_e0 << std::endl;
        std::cout << "e1   " << mainSys.e1 << std::endl;
        std::cout << "ex   " << mainSys.ex << std::endl;
        std::cout << "ey   " << mainSys.ey << std::endl;
        std::cout << "X   " << 1.75*pars.initialStretch << std::endl;
        std::cout << "Y   " << pars.initialStretch << std::endl;

        
        std::cout << "pressure   " << mainSys.P2  <<  std::endl;

        gd_solver(baseData,pars,timeStep, mainSys, false);

        if (mainSys.P2<pars.targetPressure) break;
        
        mainSys.affine_compression_triWalls(baseData, pars, - pars.deformationRate * pars.dt ,1.75*pars.initialStretch, 1*pars.initialStretch);
    }
//    long t = timeStep;
//    while(timeStep<=t+10){
//
//        std::cout << timeStep << std::endl;
//        std::cout << "***** holding *******   " << std::endl;
//        std::cout << "phi   " << mainSys.phi << "  target is  " << pars.targetPhi << std::endl;
//        std::cout << "e0   " << mainSys.e0 << "  target is  " << target_e0 << std::endl;
//        std::cout << "pressure   " << mainSys.P2  <<  std::endl;
//
//        gd_solver(baseData,pars,timeStep, mainSys, false);
//
//    }
    
    mainSys.dump_global_data(pars, timeStep, "append", "final");
    mainSys.dump_per_node(baseData, pars, timeStep);
    mainSys.dump_per_ele(baseData, pars,timeStep);
    if (pars.dumpPeriodicImagesXY){
        mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
    }
    
    timeStep=0;
    mainSys.dump_global_data(pars, timeStep, "write", "running");
    
    while(mainSys.e1 <= pars.shearTo){
        
        std::cout << timeStep << std::endl;
        std::cout << "phi   " << mainSys.phi << pars.targetPhi << std::endl;
        std::cout << "e0   " << mainSys.e0 << std::endl;
        std::cout << "e1   " << mainSys.e1 << std::endl;
        std::cout << "pressure   " << mainSys.P2  <<  std::endl;

        gd_solver(baseData,pars,timeStep, mainSys, true);

        mainSys.affine_axial_shearing_triWalls(baseData, pars, pars.deformationRate * pars.dt,1.75*pars.initialStretch, 1*pars.initialStretch);
    }
}
        
        
        
