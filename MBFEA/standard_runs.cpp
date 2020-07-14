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



void compress(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys){
    
    double refPhi = pars.Ap/(baseData.lxRef*baseData.lyRef);
    double target_e0 = - log(sqrt(pars.Ap/(pars.targetPhi*baseData.lxRef*baseData.lyRef))); // my conviension is positive e0 for compression

    
    if (pars.solver=="FIRE2"){
        
        assert(pars.FIRE_startingStrainStep <= pars.FIRE_numStrainSteps);
        
        for (long strainStep = pars.FIRE_startingStrainStep ; strainStep<= pars.FIRE_numStrainSteps ; strainStep++) {
            
            double stepPhi = refPhi + (pars.targetPhi - refPhi) * float(strainStep)/float(pars.FIRE_numStrainSteps); // the required phi of this step as a fraction of the required target phi
            double relativeStrain = - log(sqrt(pars.Ap/(stepPhi*baseData.lxRef*baseData.lyRef))) - mainSys.e0;  // strain between current and new configuration
            
            std::cout << timeStep << std::endl;
            std::cout << "phi  " << mainSys.phi << "  target is  " << pars.targetPhi << std::endl;
            std::cout << "e0  " << mainSys.e0 << "  target is  " << target_e0 << std::endl;
            
            mainSys.affine_compression(baseData, pars, relativeStrain, mainSys.ctrX, mainSys.ctrY, timeStep);
            
            fire2_solver(baseData, pars, timeStep , mainSys, strainStep);
            
            if (mainSys.maxR > pars.FIRE_RTolerance){
                 std::cout << " Failed to coverge !" << std::endl;

             }else{
                 mainSys.dump_global_data(pars, timeStep, "data", "append", "final");
                 mainSys.dump_per_node(baseData, pars, strainStep);
                 mainSys.dump_per_ele(baseData, pars,strainStep);
                 if (pars.dumpPeriodicImagesXY){
                     mainSys.dump_per_node_periodic_images_on(baseData, pars, strainStep);
                 }
             }
        }
             
            
    }else if (pars.solver=="GD") {
        
        while(1){
            
            std::cout << timeStep << std::endl;
                std::cout << "phi   " << mainSys.phi << "  target is  " << pars.targetPhi << std::endl;
                std::cout << "e0   " << mainSys.e0 << "  target is  " << target_e0 << std::endl;
                std::cout << "e1   " << mainSys.e1 << std::endl;

            gd_solver(baseData,pars,timeStep, "data", mainSys, true);
          
            if (mainSys.phi <= pars.targetPhi)
            {
                mainSys.affine_compression(baseData, pars, pars.deformationRate * pars.dt, mainSys.ctrX, mainSys.ctrY, timeStep);
            }
        }
    }
}
    

    


void stepshear(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys){
    
    mainSys.dump_global_data(pars, timeStep,"final-data", "write", "final");

    mainSys.dump_global_data(pars, timeStep,"holding-data", "write", "running");
    
    
    //hold
    while(mainSys.maxR > pars.maxForceTol){
        
         if(timeStep%pars.writeToConsoleEvery==0){
            std::cout << "**** holding before shearing **** " << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
            std::cout << timeStep << std::endl;
            std::cout << "phi   " << mainSys.phi << pars.targetPhi << std::endl;
            std::cout << "e0   " << mainSys.e0 << std::endl;
            std::cout << "e1   " << mainSys.e1 << std::endl;
            std::cout << "pressure   " << mainSys.P2  <<  std::endl;
            std::cout << "interactions  " << mainSys.nodeIinteractions + mainSys.segmentIinteractions << std::endl;

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
    
    timeStep=0;
    
    
    mainSys.affine_axial_shearing(baseData, pars, pars.targetShear ,mainSys.ctrX, mainSys.ctrY, timeStep);
    
    do {
        
         if(timeStep%pars.writeToConsoleEvery==0){
            std::cout << "**** holding after shearing **** " << "\t dt \t" << pars.dt << " \t " << pars.boundaryType << " \t " << pars.contactMethod << std::endl;
            
            std::cout << timeStep << std::endl;
            std::cout << "phi   " << mainSys.phi << pars.targetPhi << std::endl;
            std::cout << "e0   " << mainSys.e0 << std::endl;
            std::cout << "e1   " << mainSys.e1 << std::endl;
            std::cout << "pressure   " << mainSys.P2  <<  std::endl;
            std::cout << "interactions  " << mainSys.nodeIinteractions + mainSys.segmentIinteractions << std::endl;

         }
        
        gd_solver(baseData,pars,timeStep, "holding-data", mainSys, false);

        
    } while(mainSys.maxR > pars.maxForceTol);
    
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
    
    
}











////////////////////////////////////////Old way
//void stepshear(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys){
//
//   long stage=0;
//   mainSys.dump_global_data(pars, timeStep, "data","write", "final"); //open file and write cols names
//   while(stage<2){
//       if (pars.solver=="FIRE2"){
//
//           mainSys.affine_axial_shearing(baseData, pars, pars.targetShear, mainSys.ctrX, mainSys.ctrY, timeStep);
//
//           fire2_solver(baseData, pars, timeStep , mainSys, stage);
//
//            if (mainSys.maxR > pars.FIRE_RTolerance){
//                 std::cout << " Failed to coverge !" << std::endl;
//                 exit(1);
//             }else{
//                 mainSys.dump_global_data(pars, timeStep, "data", "append", "final");
//                 mainSys.dump_per_node(baseData, pars, stage);
//                 mainSys.dump_per_ele(baseData, pars,stage);
//                 if (pars.dumpPeriodicImagesXY){
//                     mainSys.dump_per_node_periodic_images_on(baseData, pars, stage);
//                 }
//                 if (pars.identifyAndDumbFacets) {
//                     if (pars.identifyAndDumbFacets) {
//                         if (pars.contactMethod=="nts"){
//                             mainSys.dump_facets(baseData, pars, timeStep);
//                         }else{
//                             mainSys.dump_facets_ntn(baseData, pars, timeStep);
//                         }
//                     }
//                 }
//
//                 stage++;
//            }
//       }else if (pars.solver=="GD") {
//
//           if(timeStep%pars.writeToConsoleEvery==0){
//           std::cout << timeStep << std::endl;
//           std::cout << "stage  " << stage << std::endl;
//           }
//           gd_solver(baseData,pars,timeStep, "data", mainSys, false);
//
//           if (mainSys.e1 < pars.targetShear) //Here tergetShear is initial nonzerovalue for shear
//           {
//               mainSys.affine_axial_shearing(baseData, pars, pars.deformationRate * pars.dt, mainSys.ctrX, mainSys.ctrY, timeStep);
//
//           }else if (mainSys.maxR <= pars.maxForceTol){
//
//               if (stage==0){
//                   mainSys.dump_global_data(pars, timeStep,"data", "write", "final");
//                   mainSys.dump_global_data(pars, timeStep,"data", "append", "final");
//                   mainSys.dump_per_node(baseData, pars, timeStep);
//                   mainSys.dump_per_ele(baseData, pars,timeStep);
//                   if (pars.dumpPeriodicImagesXY){
//                       mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
//                   }
//                   if (pars.identifyAndDumbFacets) {
//                       if (pars.identifyAndDumbFacets) {
//                           if (pars.contactMethod=="nts"){
//                               mainSys.dump_facets(baseData, pars, timeStep);
//                           }else{
//                               mainSys.dump_facets_ntn(baseData, pars, timeStep);
//                           }
//                       }
//                   }
//                   mainSys.affine_axial_shearing(baseData, pars, pars.deformationRate * pars.dt, mainSys.ctrX, mainSys.ctrY, timeStep);
//                   stage++;
//               } else if (stage==1){
//                   mainSys.dump_global_data(pars, timeStep, "data", "append", "final");
//                   mainSys.dump_per_node(baseData, pars, timeStep);
//                   mainSys.dump_per_ele(baseData, pars,timeStep);
//                   if (pars.dumpPeriodicImagesXY){
//                       mainSys.dump_per_node_periodic_images_on(baseData, pars, timeStep);
//                   }
//                   if (pars.identifyAndDumbFacets) {
//                       if (pars.identifyAndDumbFacets) {
//                           if (pars.contactMethod=="nts"){
//                               mainSys.dump_facets(baseData, pars, timeStep);
//                           }else{
//                               mainSys.dump_facets_ntn(baseData, pars, timeStep);
//                           }
//                       }
//                   }
//
//                   stage++;
//               }
//           }
//       }
//    }
//}

