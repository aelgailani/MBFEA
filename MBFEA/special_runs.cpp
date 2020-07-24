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
        F.nx = nx;
        F.ny = ny;
        F.r = r;
        F.f = f;
    }
    return F;
}


struct PowerlawRepulsion walldiscretePL(double x, double x0, double x1, double y, double y0,double y1, double sigma, double epsilon, double rcut, int numwallnodes){
    
    double dx01 = x1 - x0;
    double dy01 = y1 - y0;
     
    struct PowerlawRepulsion F;
    
     for (double s = 0; s < 1.0; s+=1.0/numwallnodes){
         
         double xs = x0 + s * dx01;
         double ys = y0 + s * dy01;
         
         double dxij = xs - x;
         double dyij = ys - y;
         double drijSq = pow(dxij,2)+ pow(dyij,2);
         double drij = sqrt(drijSq);
         
         if (drij > rcut) continue;
         
         double forceij=0.5*epsilon/sigma*12*pow((sigma/drij),13);
         double forceXij = - forceij*dxij/drij;
         double forceYij = - forceij*dyij/drij;
        
         F.fx += forceXij;
         F.fy += forceYij;
         F.energy += 0.5*epsilon*pow((sigma/drij),12);
     }
    return F;
}

void shear_special_FIRE(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys){
    float shearStep=pars.deformationRate;
    long numSteps=1;
    long step = 0;
    double refPhi = pars.Ap/(0.5*baseData.lxRef*baseData.lyRef);
    double target_e0 = - log(sqrt(refPhi/(pars.targetPhi)));
    double stepPhi = refPhi;
    
    mainSys.dump_global_data(pars, timeStep,"final-data", "write", "final");
    // initial compression
    for (long strainStep = pars.FIRE_startingStrainStep ; strainStep<= pars.FIRE_numStrainSteps ; strainStep++) {
        double prevStepPhi = stepPhi;
        stepPhi = refPhi + (pars.targetPhi - refPhi) * float(strainStep)/float(pars.FIRE_numStrainSteps); // the required phi of this step as a fraction of the required target phi
        double relativeStrain = log(sqrt(prevStepPhi/stepPhi));  // strain between current and new configuration
        
         if(timeStep%pars.writeToConsoleEvery==0){
            std::cout << "**** compressing ****" << std::endl;
            std::cout << timeStep << std::endl;
            std::cout << "phi  " << mainSys.phi * 2 << "  target is  " << pars.targetPhi << std::endl;
            std::cout << "e0  " << mainSys.e0 << "  target is  " << target_e0 << std::endl;
         }
        
        mainSys.affine_compression_triWalls(baseData, pars, relativeStrain,mainSys.ctrX, mainSys.ctrY, timeStep);
        
        fire2_solver(baseData, pars, timeStep , mainSys, strainStep);
        
        if (mainSys.L2NormResidual > pars.FIRE_RTolerance){
             std::cout << " Failed to coverge compression !" << std::endl;

         }
            else{
            mainSys.dump_per_node(baseData, pars, strainStep);
              mainSys.dump_per_ele(baseData, pars,strainStep);
              if (pars.dumpPeriodicImagesXY){
                  mainSys.dump_per_node_periodic_images_on(baseData, pars, strainStep);
              }
             if(pars.dumpSmoothenCurves){
                  mainSys.dump_smoothcurves(baseData,pars,strainStep);
              }

         }
        if (mainSys.P2<pars.targetPressure){
            
            mainSys.dump_per_node(baseData, pars, strainStep);
             mainSys.dump_per_ele(baseData, pars,strainStep);
             if (pars.dumpPeriodicImagesXY){
                 mainSys.dump_per_node_periodic_images_on(baseData, pars, strainStep);
             }
            if(pars.dumpSmoothenCurves){
                 mainSys.dump_smoothcurves(baseData,pars,strainStep);
             }
            break;
        }
    }
    
    mainSys.dump_global_data(pars, timeStep,"final-data", "append", "final");
    
    while(step<=numSteps && mainSys.e1 <= pars.shearTo){
        
         if(timeStep%pars.writeToConsoleEvery==0){
            std::cout << "**** shearing ****" << std::endl;
            std::cout << timeStep << std::endl;
            std::cout << "phi  " << mainSys.phi *2<< std::endl;
            std::cout << "e1  " << mainSys.e1 << "  target is  " << target_e0 << std::endl;
         }
        mainSys.affine_axial_shearing_triWalls(baseData, pars, shearStep,mainSys.ctrX, mainSys.ctrY, timeStep);

        
        fire2_solver(baseData, pars, timeStep , mainSys, step);
        
        if (mainSys.L2NormResidual > pars.FIRE_RTolerance){
             std::cout << " Failed to coverge !" << std::endl;
             exit(1);
//         }else{
//             mainSys.dump_global_data(pars, timeStep,"final-data","append", "final");
//             mainSys.dump_per_node(baseData, pars, step);
//             mainSys.dump_per_ele(baseData, pars,step);
//             if (pars.dumpPeriodicImagesXY){
//                 mainSys.dump_per_node_periodic_images_on(baseData, pars, step);
//             }
//             if (pars.identifyAndDumbFacets) {
//                 if (pars.contactMethod=="nts"){
//                     mainSys.dump_facets(baseData, pars, timeStep);
//                 }else{
//                     mainSys.dump_facets_ntn(baseData, pars, timeStep);
//                 }
//             }
//
//
        }
        
        step++;
    }
    if (mainSys.L2NormResidual > pars.FIRE_RTolerance){
                 std::cout << " Failed to coverge !" << std::endl;
                 exit(1);
             }else{
                 mainSys.dump_global_data(pars, timeStep,"final-data","append", "final");
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
                 }if(pars.dumpSmoothenCurves){
                     mainSys.dump_smoothcurves(baseData,pars,step);
                 }
    
    
            }
}

void shear_special_GD(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys){

    double refPhi = pars.Ap/(0.5*baseData.lxRef*baseData.lyRef);
    double target_e0 = - log(sqrt(refPhi/(pars.targetPhi)));
    
    mainSys.dump_global_data(pars, timeStep,"compression-data", "write", "running");
    // initial compression
    gd_solver(baseData,pars,timeStep,"compression-data", mainSys, true);

    
    while(mainSys.phi <= pars.targetPhi){
        
        mainSys.affine_compression_triWalls(baseData, pars, - pars.deformationRate * pars.dt ,mainSys.ctrX, mainSys.ctrY, timeStep);
         
        if(timeStep%pars.writeToConsoleEvery==0){
            std::cout << timeStep << std::endl;
            std::cout << "ctr X   " << mainSys.ctrX << std::endl;
            std::cout << "ctr Y   " << mainSys.ctrY << std::endl;
            std::cout << "height   " << mainSys.height << std::endl;
            std::cout << "base   " << mainSys.base << std::endl;
            std::cout << "phi   " << mainSys.phi << "  target is  " << pars.targetPhi << std::endl;
            std::cout << "e0   " << mainSys.e0 << "  target is  " << target_e0 << std::endl;
            std::cout << "e1   " << mainSys.e1 << std::endl;
            std::cout << "pressure   " << mainSys.P2  <<  std::endl;
            std::cout << "interactions  " << mainSys.nodeIinteractions + mainSys.segmentIinteractions << std::endl;
            std::cout << "maxForce  " << mainSys.maxR << std::endl;
            std::cout << "meanForce  " << mainSys.avgR << std::endl;
            std::cout << "L2NormR  " << mainSys.L2NormResidual << std::endl;
            std::cout << "\n" << std::endl;


         }
        
        gd_solver(baseData,pars,timeStep,"compression-data", mainSys, true);

        if (mainSys.P2<pars.targetPressure) break;
        
       
    }

    mainSys.dump_global_data(pars, timeStep, "compression-data", "append", "final");
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
    timeStep=1;
    
    mainSys.dump_global_data(pars, timeStep, "shearing-data","write", "running");
    
    while(mainSys.e1 <= pars.shearTo){
        
        mainSys.affine_axial_shearing_triWalls(baseData, pars, pars.deformationRate * pars.dt,mainSys.ctrX, mainSys.ctrY, timeStep);
         
        if(timeStep%pars.writeToConsoleEvery==0){
            std::cout << timeStep << std::endl;
            std::cout << "ctr X   " << mainSys.ctrX << std::endl;
            std::cout << "ctr Y   " << mainSys.ctrY << std::endl;
            std::cout << "height   " << mainSys.height << std::endl;
            std::cout << "base   " << mainSys.base << std::endl;
            std::cout << "phi   " << mainSys.phi << pars.targetPhi << std::endl;
            std::cout << "e0   " << mainSys.e0 << std::endl;
            std::cout << "e1   " << mainSys.e1 << std::endl;
            std::cout << "pressure   " << mainSys.P2  <<  std::endl;
            std::cout << "interactions  " << mainSys.nodeIinteractions + mainSys.segmentIinteractions << std::endl;
            std::cout << "meanForce  " << mainSys.avgR << std::endl;
            std::cout << "L2NormR  " << mainSys.L2NormResidual << std::endl;
            std::cout << "\n" << std::endl;


         }
        
        gd_solver(baseData,pars,timeStep,"shearing-data", mainSys, true);

        
    }
}
        
        
        
void shear_special_stepGD(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys){

    double refPhi = pars.Ap/(0.5*baseData.lxRef*baseData.lyRef);
    double target_e0 = - log(sqrt(refPhi/(pars.targetPhi)));
    mainSys.dump_global_data(pars, timeStep,"compression-data", "write", "running");
    mainSys.dump_global_data(pars, timeStep,"final-data", "write", "final");
    // initial compression
    gd_solver(baseData,pars,timeStep,"compression-data", mainSys, true);

    while(mainSys.phi <= pars.targetPhi){
        
        mainSys.affine_compression_triWalls(baseData, pars, - pars.deformationRate * pars.dt ,mainSys.ctrX, mainSys.ctrY, timeStep);
        
         if(timeStep%pars.writeToConsoleEvery==0){
            std::cout << timeStep << std::endl;
            std::cout << "ctr X   " << mainSys.ctrX << std::endl;
            std::cout << "ctr Y   " << mainSys.ctrY << std::endl;
            std::cout << "height   " << mainSys.height << std::endl;
            std::cout << "base   " << mainSys.base << std::endl;
            std::cout << "phi   " << mainSys.phi << "  target is  " << pars.targetPhi << std::endl;
            std::cout << "e0   " << mainSys.e0 << "  target is  " << target_e0 << std::endl;
            std::cout << "e1   " << mainSys.e1 << std::endl;
            std::cout << "pressure   " << mainSys.P2  <<  std::endl;
            std::cout << "interactions  " << mainSys.nodeIinteractions + mainSys.segmentIinteractions << std::endl;
            std::cout << "meanForce  " << mainSys.avgR << std::endl;
            std::cout << "L2NormR  " << mainSys.L2NormResidual << std::endl;
            std::cout << "\n" << std::endl;

         }
        
        gd_solver(baseData,pars,timeStep, "compression-data",mainSys, false);

        if (mainSys.P2<pars.targetPressure) break;
        
       
    }

    mainSys.dump_global_data(pars, timeStep, "final-data", "append", "final");
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
    timeStep=0;
    
    mainSys.dump_global_data(pars, timeStep,"holding-data", "write", "running");
    
    
    //hold
    while(mainSys.L2NormResidual > pars.maxForceTol){
        
         if(timeStep%pars.writeToConsoleEvery==0){
            std::cout << "**** holding before shearing **** " << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
            std::cout << timeStep << std::endl;
            std::cout << "ctr X   " << mainSys.ctrX << std::endl;
            std::cout << "ctr Y   " << mainSys.ctrY << std::endl;
            std::cout << "height   " << mainSys.height << std::endl;
            std::cout << "base   " << mainSys.base << std::endl;
            std::cout << "phi   " << mainSys.phi << pars.targetPhi << std::endl;
            std::cout << "e0   " << mainSys.e0 << std::endl;
            std::cout << "e1   " << mainSys.e1 << std::endl;
            std::cout << "pressure   " << mainSys.P2  <<  std::endl;
            std::cout << "interactions  " << mainSys.nodeIinteractions + mainSys.segmentIinteractions << std::endl;
            std::cout << "meanForce  " << mainSys.avgR << std::endl;
            std::cout << "L2NormR  " << mainSys.L2NormResidual << std::endl;
            std::cout << "\n" << std::endl;


            
         }
        
        gd_solver(baseData,pars,timeStep, "holding-data", mainSys, false);

        
    }
    
    
    mainSys.dump_global_data(pars, timeStep, "final-data","append", "final");
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
    timeStep=0;
    
    
    mainSys.affine_axial_shearing_triWalls(baseData, pars, pars.shearTo ,mainSys.ctrX, mainSys.ctrY, timeStep);
    
    do {
        
         if(timeStep%pars.writeToConsoleEvery==0){
            std::cout << "**** holding after shearing **** " << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
            
            std::cout << timeStep << std::endl;
            std::cout << "ctr X   " << mainSys.ctrX << std::endl;
            std::cout << "ctr Y   " << mainSys.ctrY << std::endl;
            std::cout << "height   " << mainSys.height << std::endl;
            std::cout << "base   " << mainSys.base << std::endl;
            std::cout << "phi   " << mainSys.phi << pars.targetPhi << std::endl;
            std::cout << "e0   " << mainSys.e0 << std::endl;
            std::cout << "e1   " << mainSys.e1 << std::endl;
            std::cout << "pressure   " << mainSys.P2  <<  std::endl;
            std::cout << "interactions  " << mainSys.nodeIinteractions + mainSys.segmentIinteractions << std::endl;
            std::cout << "meanForce  " << mainSys.avgR << std::endl;
            std::cout << "L2NormR  " << mainSys.L2NormResidual << std::endl;
            std::cout << "\n" << std::endl;


         }
        
        gd_solver(baseData,pars,timeStep, "holding-data", mainSys, false);

        
    } while(mainSys.L2NormResidual > pars.maxForceTol);
    
        mainSys.dump_global_data(pars, timeStep, "final-data","append", "final");
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
