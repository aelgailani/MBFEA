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
#include "Parameters.hpp"
#include "Configuration.hpp"
#include "BaseSysData.hpp"

void Configuration::contact_forces(const BaseSysData& baseData, const Parameters& pars, bool Hessian, const long& timeStep)
{
    maxInterference = 0;
    segmentIinteractions = 0;
    nodeIinteractions = 0;
    
    // remember imagesMargies = 0 if walls are used 
    numXCells = floor(lxNew*(1+2*pars.imagesMargin)/pars.verletCellCutoff);
    numYCells = floor(lyNew*(1+2*pars.imagesMargin)/pars.verletCellCutoff);
    verletCellSizeX  = lxNew*(1+2*pars.imagesMargin)/numXCells;
    verletCellSizeY = lyNew*(1+2*pars.imagesMargin)/numYCells;
    
    
    nodesLinkedList.resize(numXCells*numYCells+baseData.numSurfaceNodes,-1);
    
    //clear factes map only if required to not hinder the performance; since factes are onlny for post processing
    if (timeStep % pars.dumpEvery == 0 && pars.identifyAndDumbFacets) {
        facets.clear();
        facets_ntn.clear();
    }
    
    if (pars.contactMethod=="nts"){
        
        // clear previous data
        surNodes_mMesh1.resize(baseData.numSurfaceNodes,-1);
        surNodes_mMesh2.resize(baseData.numSurfaceNodes,-1);
        surNodes_mMesh3.resize(baseData.numSurfaceNodes,-1);
        surNodes_mSegment1.resize(baseData.numSurfaceNodes,-1);
        surNodes_mSegment2.resize(baseData.numSurfaceNodes,-1);
        surNodes_mSegment3.resize(baseData.numSurfaceNodes,-1);
        surNodes_mPart1.resize(baseData.numSurfaceNodes,-1);
        surNodes_mPart2.resize(baseData.numSurfaceNodes,-1);
        surNodes_mPart3.resize(baseData.numSurfaceNodes,-1);
        surNodes_gap1.resize(baseData.numSurfaceNodes,-1);
        surNodes_gap2.resize(baseData.numSurfaceNodes,-1);
        surNodes_gap3.resize(baseData.numSurfaceNodes,-1);
        
        
        if (pars.segmentCellMethod ==1){
            
            if (not pars.reversibleMasterSlaveRole){
                slaveMaster.resize(baseData.numMeshes,baseData.numMeshes);
                slaveMaster.fill(-1);
            }
            
            segmentsLinkedList_1.resize(numXCells*numYCells+baseData.numSurfaceNodes,-1);
            update_cells_1(baseData, pars);
            detect_nts_contacts_single_point_method(baseData,pars);
            
        }else if (pars.segmentCellMethod ==2){
            
            if (not pars.reversibleMasterSlaveRole){
                slaveMaster.resize(baseData.numMeshes,baseData.numMeshes);
                slaveMaster.fill(-1);
            }
            segmentsLinkedList_2.resize(baseData.numSurfaceNodes, 5);
            segmentsLinkedList_2.fill(-2);
            cellsHeads.resize(numXCells*numYCells, 2);
            cellsHeads.fill(-1);
            update_cells_2(baseData, pars);
            detect_nts_contacts_two_points_method(baseData,pars);
            
        }else{
            std::cout << "Please specify segmentCellMethod. Either 1 or 2. " << std::endl;;
            exit(1);
        }
        if (pars.ntsPenaltyMethod=="harmonic"){
            apply_nts_harmonic_penalty(baseData, pars, surNodes_mMesh1, surNodes_mSegment1, surNodes_mPart1, surNodes_gap1,Hessian, timeStep);
            apply_nts_harmonic_penalty(baseData, pars, surNodes_mMesh2, surNodes_mSegment2, surNodes_mPart2, surNodes_gap2, Hessian, timeStep);
            apply_nts_harmonic_penalty(baseData, pars, surNodes_mMesh3, surNodes_mSegment3, surNodes_mPart1, surNodes_gap3, Hessian, timeStep);
        }else if (pars.ntsPenaltyMethod=="powerlaw"){
           apply_nts_powerlaw_penalty(baseData, pars, surNodes_mMesh1, surNodes_mSegment1, surNodes_mPart1, surNodes_gap1,Hessian, timeStep);
            apply_nts_powerlaw_penalty(baseData, pars, surNodes_mMesh2, surNodes_mSegment2, surNodes_mPart2, surNodes_gap2, Hessian, timeStep);
            apply_nts_powerlaw_penalty(baseData, pars, surNodes_mMesh3, surNodes_mSegment3, surNodes_mPart1, surNodes_gap3, Hessian, timeStep);
            
        }else{
            std::cout << "Please specify a valid ntsPenaltyMethod. Either powerlaw or harmonic . " << std::endl;;
            exit(1);
        }
        
    }else if (pars.contactMethod=="ntn"){
       
        update_cells_1(baseData, pars);
        apply_ntn_repulsions( baseData, pars, Hessian, timeStep);
        
    }else if (pars.contactMethod=="gntn"){
       
        update_cells_1(baseData, pars);
        apply_ghost_nodes_to_node_repulsions( baseData, pars, Hessian, timeStep);

        
    }else if (pars.contactMethod=="gntn2"){
       
        if (not pars.reversibleMasterSlaveRole){
            slaveMaster.resize(baseData.numMeshes,baseData.numMeshes);
            slaveMaster.fill(-1);
        }
        segmentsLinkedList_2.resize(baseData.numSurfaceNodes, 5);
        segmentsLinkedList_2.fill(-2);
        cellsHeads.resize(numXCells*numYCells, 2);
        cellsHeads.fill(-1);
        update_cells_2(baseData, pars);
        
        apply_ghost_nodes_to_node_repulsions_v2( baseData, pars, Hessian, timeStep);

        
    }else{
        std::cout << "Please specify a valid contactMethod. Either ntn or nts . " << std::endl;;
        exit(1);
        
    }
    
    

    
    
    
}
