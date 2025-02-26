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
#include "utility_functions.hpp"


void Configuration::dump_facets(const BaseSysData& baseData, const Parameters& pars, long& timeStep){
    
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    if (pars.runMode == "stepShear" || pars.runMode == "surfaceShear" || pars.runMode == "continuousShear"){
        if (pars.outputSubFolder=="none"){
            myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingStepNum))+"/facets-"+step+".txt").c_str());
        }else{
            myfile.open ((pars.outputFolderName +"/"+ pars.outputSubFolder+"/facets-"+step+".txt").c_str());
        }
        

    }else{
        myfile.open (pars.outputFolderName+"/facets-"+step+".txt");
    }
    myfile << "Basic_data:" << std::endl;
    myfile << "timeStep" << "\t" << timeStep << std::endl;
    myfile << "\n" << std::endl;
    myfile << "Facets_data: (in the smNodes, the first entery is a slave node followed by its masters or master nodes or node, then the next slave node and so forth. Notice that slave nodes never duplicate in a row but master nodes can be shared and hence duplicated. Each slave node has two master nodes at max at one master mesh and a minimum of one.)" << std::endl;
    myfile
    <<"sMesh" << std::setw(7)
    << "mMesh" << std::setw(7)
    << "smNodes" <<  std::endl;

    for(auto& key : facets) {
        myfile << std::to_string(key.first.first) << std::setw(7) << std::to_string(key.first.second) << std::setw(7);
        
        for (auto& node: key.second) {
            myfile << std::to_string(node) << std::setw(7);
        }
        myfile << std::setw(-7)<< std::endl;
    }
    myfile << "EOF";
    myfile.close();

}

void Configuration::dump_facets_ntn(const BaseSysData& baseData, const Parameters& pars, long& timeStep){
    
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    if (pars.runMode == "stepShear" || pars.runMode == "surfaceShear"  || pars.runMode == "continuousShear"){
        if (pars.outputSubFolder=="none"){
                myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingStepNum))+"/facets-"+step+".txt").c_str());
            }else{
                myfile.open ((pars.outputFolderName +"/"+ pars.outputSubFolder+"/facets-"+step+".txt").c_str());
            }

    }else{
        myfile.open (pars.outputFolderName+"/facets-"+step+".txt");
    }
    myfile << "Basic_data:" << std::endl;
    myfile << "timeStep" << "\t" << timeStep << std::endl;
    myfile << "\n" << std::endl;
    
    
    if(pars.contactMethod == "ntn"){
        myfile << "Facets_data:" << std::endl;
        myfile
        <<"iMesh" << std::setw(20)
        << "jMesh" << std::setw(20)
        << "inode" << std::setw(20)
        << "jnode" << std::setw(20)
        << "d" << std::setw(20)
        << "f" << std::setw(20)
        << "ix" << std::setw(20)
        << "jx" << std::setw(20)
        << "iy" << std::setw(20)
        << "jy" << std::setw(20)
        << "fx" << std::setw(20)
        << "fy" <<  std::endl;
    }else if (pars.contactMethod == "gntn"){
        myfile << "Facets_data: " << std::endl;
        myfile
        << "iMesh" << std::setw(20)
        << "jMesh" << std::setw(20)
        << "inode" << std::setw(20)
        << "s" << std::setw(20)
        << "d" << std::setw(20)
        << "f" << std::setw(20)
        << "ix" << std::setw(20)
        << "sx" << std::setw(20)
        << "iy" << std::setw(20)
        << "sy" << std::setw(20)
        << "fx" << std::setw(20)
        << "fy" <<  std::setw(20)
        << " there are " << std::to_string(pars.gntn_NGhostNodes) << "  ghost nodes on each segment" << std::endl;
    }
    

    for(auto& key : facets_ntn) {
        myfile << std::to_string(key.first.first) << std::setw(20) << std::to_string(key.first.second) << std::setw(20);
        
        for (auto& node: key.second) {
            myfile << std::to_string(node) << std::setw(20);
        }
        myfile << std::setw(-20)<< std::endl;
    }
    myfile << "EOF";
    myfile.close();

}


void Configuration::dump_global_data(const Parameters& pars, const long& timeStep, std::string name, std::string mode, std::string purpose){
    int spacing =26;
    std::string fname ;
    std::string first = std::to_string(timeStep/pars.splitDataEvery*pars.splitDataEvery);
    std::string last =std::to_string(timeStep/pars.splitDataEvery*pars.splitDataEvery+pars.splitDataEvery-1);
    
    
    if (purpose=="running"){
        fname = "/"+name+"-steps-"+first+"-"+last+".txt";
    }else if (purpose=="final"){
        fname = "/final-data.txt";
    }
   
    if (mode=="write"){
        std::ofstream dataFile;
        if (pars.runMode == "stepShear" || pars.runMode == "surfaceShear" || pars.runMode == "continuousShear"){
        
        
            if (pars.outputSubFolder=="none"){
                dataFile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingStepNum))).c_str()+fname);
            }else{
                dataFile.open ((pars.outputFolderName +"/"+ pars.outputSubFolder).c_str()+fname);
            }

        }else{
            dataFile.open (pars.outputFolderName+fname);
        }
        dataFile
        <<  "step"  << std::setw(25)
        <<  "totalEnergy"  << std::setw(spacing)
        <<  "contactEnergy"  << std::setw(spacing)
        <<  "shearVirial"  << std::setw(spacing)
        <<  "pressureVirial"  << std::setw(spacing)
        <<  "maxResidual"  << std::setw(spacing)
        <<  "meanResidual"  << std::setw(spacing)
        <<  "residualL2Norm"  << std::setw(spacing)
        
        <<  "pressure_BW"  << std::setw(spacing)
        <<  "pressure_BB"  << std::setw(spacing)
        <<  "KWpressure_BW"  << std::setw(spacing)
        
        
        <<  "shearStress_BW"  << std::setw(spacing)
        <<  "shearStress_BB"  << std::setw(spacing)
        <<  "shearStress_BB_CSxy"  << std::setw(spacing)
        <<  "KWshearStress_BW"  << std::setw(spacing)
        
        <<  "boxArea"  << std::setw(spacing)
        <<  "materialArea"  << std::setw(spacing)
        <<  "phi"  << std::setw(spacing)
        <<  "phi_B"  << std::setw(spacing)
        <<  "e0_W"  << std::setw(spacing)
        <<  "e1_W" << std::setw(spacing)
        <<  "e0_B"  << std::setw(spacing)
        <<  "e1_B" << std::setw(spacing)
        <<  "Fh_B" << std::setw(spacing)
        <<  "Fv_B" << std::setw(spacing)
        <<  "Fh_B_CSxy" << std::setw(spacing)
        <<  "Fv_B_CSxy" << std::setw(spacing)
        <<  "LXref_B" << std::setw(spacing)
        <<  "LYref_B" << std::setw(spacing)
        <<  "LX_B" << std::setw(spacing)
        <<  "LY_B" << std::setw(spacing)
        
        
        <<  "DPOverDe0_BW" << std::setw(spacing)
        <<  "DSOverDe1_BW" << std::setw(spacing)
    
        <<  "DPOverDe0_BB" << std::setw(spacing)
        <<  "DSOverDe1_BB";
        if (pars.boundaryType=="walls"){
            dataFile
            << std::setw(spacing)
            <<  "pressure_WW"  << std::setw(spacing)
            <<  "shearStress_WW"  << std::setw(spacing)
            <<  "Fh_W" << std::setw(spacing)
            <<  "Fv_W" << std::setw(spacing)
            
            
            <<  "topForce" << std::setw(spacing)
            <<  "botForce" << std::setw(spacing)
            <<  "rightForce" << std::setw(spacing)
            <<  "leftForce";
        }
        if (pars.solver=="GD"){
            dataFile
            << std::setw(spacing)
            <<  "dt" << std::setw(spacing)
            <<  "defRate" ;
        }
        if (pars.contactMethod=="nts"){
            dataFile
            << std::setw(spacing)
            <<  "penalty";
        }
        dataFile
        << std::setw(spacing)
        <<  "interactions" << std::endl;
        dataFile.close();
        
        
    }else if (mode=="append"){
        std::ofstream dataFile;
        if (pars.runMode == "stepShear" || pars.runMode == "surfaceShear" || pars.runMode == "continuousShear"){

            if (pars.outputSubFolder=="none"){
                dataFile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingStepNum))).c_str()+fname, std::ios_base::app);
            }else{
                dataFile.open ((pars.outputFolderName +"/"+ pars.outputSubFolder).c_str()+fname,  std::ios_base::app);
            }
        
        }else{
            dataFile.open (pars.outputFolderName+fname,  std::ios_base::app);
        }
        dataFile << std::setprecision(15)
        <<  timeStep  << std::setw(30)
        <<  totalEnergy  << std::setw(spacing)
        <<  contactsEnergy  << std::setw(spacing)
        <<  shearVirial  << std::setw(spacing)
        <<  pressureVirial  << std::setw(spacing)
        <<  maxR  << std::setw(spacing)
        <<  avgR  << std::setw(spacing)
        <<  L2NormResidual  << std::setw(spacing)
        
        <<  P_BW  << std::setw(spacing)
        <<  P_BB  << std::setw(spacing)
        <<  (KWoodYY+KWoodXX)/(2*LX_W*LY_W)  << std::setw(spacing)
        <<  S_BW  << std::setw(spacing)
        <<  S_BB  << std::setw(spacing)
        <<  S_BB_CSxy  << std::setw(spacing)
        <<  (KWoodYY-KWoodXX)/(2*LX_W*LY_W)   << std::setw(spacing)
        <<  A_W  << std::setw(spacing)
        <<  A_B  << std::setw(spacing)
        <<  phi  << std::setw(spacing)
        <<  phi_B  << std::setw(spacing)
        <<  e0_W  << std::setw(spacing)
        <<  e1_W << std::setw(spacing)
        <<  e0_B  << std::setw(spacing)
        <<  e1_B << std::setw(spacing)
        <<  Fh_B << std::setw(spacing)
        <<  Fv_B << std::setw(spacing)
        <<  Fh_B_CSxy << std::setw(spacing)
        <<  Fv_B_CSxy << std::setw(spacing)
        <<  refLX_B << std::setw(spacing)
        <<  refLY_B << std::setw(spacing)
        <<  LX_B << std::setw(spacing)
        <<  LY_B << std::setw(spacing)
        
        <<  DPOverDe0_BW << std::setw(spacing)
        <<  DSOverDe1_BW << std::setw(spacing)
        <<  DPOverDe0_BB << std::setw(spacing)
        <<  DSOverDe1_BB ;
        if (pars.boundaryType=="walls"){
            dataFile
            << std::setw(spacing)
            <<  P_WW << std::setw(spacing)
            <<  S_WW  << std::setw(spacing)
            <<  Fh_W << std::setw(spacing)
            <<  Fv_W << std::setw(spacing)
            <<  wallForceTop.sum() << std::setw(spacing)
            <<  wallForceBottom.sum() << std::setw(spacing)
            <<  wallForceRight.sum() << std::setw(spacing)
            <<  wallForceLeft.sum();
        }
        if (pars.solver=="GD"){
            dataFile
            << std::setw(spacing)
            <<  pars.dt << std::setw(spacing)
            <<  pars.deformationRate;
        }
        if (pars.contactMethod=="nts"){
            dataFile
            << std::setw(spacing)
            <<  pars.ntsHarmonicPenaltyStiffness;
        }
        dataFile
        << std::setw(spacing)
        <<  nodeIinteractions << std::endl;
        dataFile.close();
        
        
    }
    
   
}


void Configuration::dump_per_node(const BaseSysData& baseData, const Parameters& pars, long& timeStep){
    int spacing =26;
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    if (pars.runMode == "stepShear" || pars.runMode == "surfaceShear"  || pars.runMode == "continuousShear"){
    
        if (pars.outputSubFolder=="none"){
                myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingStepNum))+"/data-per-node-"+step+".txt").c_str());
            }else{
                myfile.open ((pars.outputFolderName +"/"+ pars.outputSubFolder+"/data-per-node-"+step+".txt").c_str());
            }
    
    }else{
        myfile.open (pars.outputFolderName+"/data-per-node-"+step+".txt");
    }
    
    myfile << "Basic_data:" << std::endl;
    myfile << "step" << "\t" << timeStep << std::endl;
    myfile << "timeStep" << "\t" << timeStep << std::endl;
    myfile << "number_of_nodes" << "\t" << baseData.numOriginalNodes << std::endl;
    myfile << "kTOverOmega" << "\t" << pars.kTOverOmega << std::endl;
    myfile << "Ap" << "\t" << pars.Ap <<  std::endl;
    myfile << "phi" << "\t" << phi << std::endl;
    myfile << "phi_B" << "\t" << phi_B << std::endl;
    myfile << "pressure_BW" << "\t" << P_BW << std::endl;
    myfile << "pressure_BB" << "\t" << P_BB << std::endl;
    myfile << "shear_stress_BW" << "\t" << S_BW << std::endl;
    myfile << "shear_stress_BB" << "\t" << S_BB << std::endl;
    myfile << "Fh_B" << "\t" << Fh_B << std::endl;
    myfile << "Fv_B" << "\t" << Fv_B << std::endl;
    myfile << "Fh_B_CSxy" << "\t" << Fh_B_CSxy << std::endl;
    myfile << "Fv_B_CSxy" << "\t" << Fv_B_CSxy << std::endl;
    myfile << "shear_stress_BB" << "\t" << S_BB << std::endl;
    myfile << "e0_strain_W" << "\t" << e0_W << std::endl;
    myfile << "e0_strain_B" << "\t" << e0_B << std::endl;
    myfile << "e1_strain_W" << "\t" << e1_W << std::endl;
    myfile << "e1_strain_B" << "\t" << e1_B << std::endl;
    myfile << "tol_max_residual" << "\t" << pars.maxForceTol << std::endl;
    myfile << "max_residual" << "\t" << maxR << std::endl;
    myfile << "mean_residual" << "\t" << avgR << std::endl;
    myfile << "L2Residual" << "\t" << L2NormResidual << std::endl;
    myfile << "verlet_cell_cutoff" << "\t" << pars.verletCellCutoff << std::endl;
    myfile << "boundariesType" << "\t" << pars.boundaryType << std::endl;
    if (pars.boundaryType=="periodic"){
        myfile << "images_margin" << "\t" << pars.imagesMargin << std::endl;
    }
    myfile << "ref_box_LXLY_B" << "\t" << refLX_B << "\t" << refLY_B << std::endl;
    myfile << "cur_box_LXLY_B" << "\t" << LX_B << "\t" << LY_B << std::endl;
    myfile << "ref_box_LXLY_W" << "\t" << refLX_W << "\t" << refLY_W << std::endl;
    myfile << "cur_box_LXLY_W" << "\t" << LX_W << "\t" << LY_W << std::endl;
    myfile << "ref_box_LRBT" << "\t" << pars.initLeftPos << "\t" << pars.initRightPos << "\t" << pars.initBotPos << "\t" << pars.initTopPos << std::endl;
    myfile << "cur_box_LRBT" << "\t" << leftPos << "\t" << rightPos << "\t" << botPos << "\t" << topPos<<  "\n--------------------------"   << std::endl;
    
    
    myfile << "solver" << "\t" << pars.solver << std::endl;
       if (pars.solver=="GD"){
           myfile << "dt" << "\t" << pars.dt << std::endl;
           myfile << "deformation_rate" << "\t" << pars.deformationRate   << std::endl;

       }
    
    myfile << "contact_method" << "\t" << pars.contactMethod << std::endl;
    if (pars.contactMethod=="ntn" || pars.contactMethod=="gntn"){
        myfile << "sigma" << "\t" << pars.ntnRadius <<  std::endl;
        myfile << "RcutOverSigma" << "\t" << pars.ntnPLRcutoffOverRadius <<  "\n--------------------------"   <<  std::endl;
        if (pars.contactMethod=="gntn"){
            myfile << "ghost_nodes_per_segment" << "\t" << pars.gntn_NGhostNodes <<  "\n--------------------------"   <<  std::endl;
        }

    }else if (pars.contactMethod=="nts"){
        myfile << "max_penetration" << "\t" << maxInterference << std::endl;
        myfile << "penalty_stiffness" << "\t" << pars.ntsHarmonicPenaltyStiffness << std::endl;
        if (pars.smoothCorners){
            myfile << "Hermit_alpha" << "\t" << pars.alpha_HermitPol <<  "\n--------------------------"   <<  std::endl;
        }
       
    }


    
    myfile << "Nodes_data:" << std::endl;
    myfile
    << "id"  << std::setw(spacing)
    <<  "x"  << std::setw(spacing)
    <<  "y"  << std::setw(spacing)
    <<  "fx"  << std::setw(spacing)
    <<  "fy"  << std::setw(spacing)
    <<  "homVx"  << std::setw(spacing)
    <<  "homVy"  << std::setw(spacing)
    <<  "totVx"  << std::setw(spacing)
    <<  "totVy"  << std::setw(spacing)
    <<  "contactFx"  << std::setw(spacing)
    <<  "contactFy";
    if (pars.boundaryType=="walls"){
        myfile
        << std::setw(spacing)
        <<  "wallForceTop"  << std::setw(spacing)
        <<  "wallForceBot"  << std::setw(spacing)
        <<  "wallForceRight"  << std::setw(spacing)
        <<  "wallForceLeft"  << std::endl;
    }else{
        myfile << std::endl;
    }
        
    
    for (int i=0; i<curPosX.size();i++ ) {
        myfile
        <<  std::setprecision(15)
        << i << std::setw(spacing)
        << curPosX[i] << std::setw(spacing)
        << curPosY[i] << std::setw(spacing)
        << forceX[i] << std::setw(spacing)
        << forceY[i] << std::setw(spacing)
        << homVx[i] << std::setw(spacing)
        << homVy[i] << std::setw(spacing)
        << homVx[i]+forceX[i] << std::setw(spacing)
        << homVy[i]+forceY[i] << std::setw(spacing)
        << surfaceForceX[i] << std::setw(spacing)
        << surfaceForceY[i];
        if (pars.boundaryType=="walls"){
            myfile
            << std::setw(spacing)
            <<  wallForceTop[i] << std::setw(spacing)
            <<  wallForceBottom[i]  << std::setw(spacing)
            <<  wallForceRight[i]  << std::setw(spacing)
            <<  wallForceLeft[i]  << std::endl;
        }else{
            myfile << std::endl;
        }
    }
    
    myfile.close();

}

void Configuration::dump_per_node_periodic_images_on(const BaseSysData& baseData, const Parameters& pars, long& timeStep){
    int spacing =26;
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    if (pars.runMode == "stepShear" || pars.runMode == "surfaceShear" || pars.runMode == "continuousShear"){
    
            if (pars.outputSubFolder=="none"){
                myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingStepNum))+"/data-per-node-periodic-"+step+".txt").c_str());
            }else{
                myfile.open ((pars.outputFolderName +"/"+ pars.outputSubFolder+"/data-per-node-periodic-"+step+".txt").c_str());
            }

    }else{
        myfile.open (pars.outputFolderName+"/data-per-node-periodic-"+step+".txt");
    }
    myfile << "Basic_data:" << std::endl;
    myfile << "timeStep" << "\t" << timeStep << std::endl;
    myfile << "number_of_augmented_nodes" << "\t" << baseData.numNodes << std::endl;
    myfile << "images_margin" << "\t" << pars.imagesMargin << std::endl;
    myfile << "effective_simulation_box_left_boundary" << "\t" << leftPos - pars.imagesMargin*newLX_W << std::endl;
    myfile << "effective_simulation_box_bot_boundary" << "\t" << botPos - pars.imagesMargin*newLY_W << std::endl;
    myfile << "numCellsX" << "\t" << numXCells << std::endl;
    myfile << "numCellsY" << "\t" << numYCells << std::endl;
    myfile << "cellSizeX" << "\t" << verletCellSizeX << std::endl;
    myfile << "cellSizeY" << "\t" << verletCellSizeY <<"\n" << std::endl;
   
    myfile << "Nodes_data:" << std::endl;
    myfile
    << "id"  << std::setw(spacing)
    <<  "x"  << std::setw(spacing)
    <<  "y"  << std::setw(spacing)
    <<  "fx"  << std::setw(spacing)
    <<  "fy"  << std::setw(spacing)
    <<  "contactFx"  << std::setw(spacing)
    <<  "contactFy"  << std::endl;
    
    for (int i=0; i<augmentedCurPosX.size();i++ ) {
        myfile
        <<  std::setprecision(15)
        << i << std::setw(spacing)
        << augmentedCurPosX[i] << std::setw(spacing)
        << augmentedCurPosY[i] << std::setw(spacing)
        << forceX[i % baseData.numOriginalNodes] << std::setw(spacing)
        << forceY[i % baseData.numOriginalNodes] << std::setw(spacing)
        << surfaceForceX[i % baseData.numOriginalNodes] << std::setw(spacing)
        << surfaceForceY[i % baseData.numOriginalNodes] << std::endl;
    }
    myfile.close();
    
}


void Configuration::dump_per_ele(const BaseSysData& baseData, const Parameters& pars, long& timeStep){
    
    int spacing =26;
    // calculated strains
    DtotVxDx = gradX * (homVx+forceX);
    DtotVxDy = gradY * (homVx+forceX);
    DtotVyDx = gradX * (homVy+forceY);
    DtotVyDy = gradY * (homVy+forceY);
    
    DhomVxDx = gradX * (homVx);
    DhomVxDy = gradY * (homVx);
    DhomVyDx = gradX * (homVy);
    DhomVyDy = gradY * (homVy);
    
    DVxDx = gradX * (forceX);
    DVxDy = gradY * (forceX);
    DVyDx = gradX * (forceY);
    DVyDy = gradY * (forceY);
    
    
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    if (pars.runMode == "stepShear" || pars.runMode == "surfaceShear"  || pars.runMode == "continuousShear"){
    

        if (pars.outputSubFolder=="none"){
                myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingStepNum))+"/data-per-ele-"+step+".txt").c_str());
            }else{
                myfile.open ((pars.outputFolderName +"/"+ pars.outputSubFolder+"/data-per-ele-"+step+".txt").c_str());
            }

    }else{
        myfile.open (pars.outputFolderName+"/data-per-ele-"+step+".txt");
    }
    myfile << "Basic_data:" << std::endl;
    myfile << "timeStep" << "\t" << timeStep << std::endl;
    myfile << "number_of_elements" << "\t" << baseData.numElements << std::endl;
    myfile << "internal_energy" << "\t" <<  std::setprecision(9)<< internalEnergy << std::endl;
    myfile << "contacts_energy" << "\t" <<  std::setprecision(9)<< contactsEnergy << std::endl;
    myfile << "total_energy" << "\t" <<  std::setprecision(9)<< totalEnergy <<"\n"<< std::endl;
    myfile << "Elements_data:" << std::endl;
    myfile
    <<"id" << std::setw(spacing)
    << "refArea" << std::setw(spacing)
    << "areaRatio" << std::setw(spacing)
    << "FXX" << std::setw(spacing)
    << "FXY" << std::setw(spacing)
    << "FYX" << std::setw(spacing)
    << "FYY" << std::setw(spacing)
    << "deltaFXX" << std::setw(spacing)
    << "deltaFXY" << std::setw(spacing)
    << "deltaFYX" << std::setw(spacing)
    << "deltaFYY" << std::setw(spacing)
    << "PK1StressXX" << std::setw(spacing)
    << "PK1StressXY" << std::setw(spacing)
    << "PK1StressYX" << std::setw(spacing)
    << "PK1StressYY" << std::setw(spacing)
    << "CStressXX" << std::setw(spacing)
    << "CStressXY" << std::setw(spacing)
    << "CStressYX" << std::setw(spacing)
    << "CStressYY" << std::setw(spacing)
    << "deltaCStressXX" << std::setw(spacing)
    << "deltaCStressXY" << std::setw(spacing)
    << "deltaCStressYX" << std::setw(spacing)
    << "deltaCStressYY" << std::setw(spacing)
    << "DtotVxDx" << std::setw(spacing)
    << "DtotVxDy" << std::setw(spacing)
    << "DtotVyDx" << std::setw(spacing)
    << "DtotVyDy" << std::setw(spacing)
    << "DhomVxDx" << std::setw(spacing)
    << "DhomVxDy" << std::setw(spacing)
    << "DhomVyDx" << std::setw(spacing)
    << "DhomVyDy" << std::setw(spacing)
    << "DVxDx" << std::setw(spacing)
    << "DVxDy" << std::setw(spacing)
    << "DVyDx" << std::setw(spacing)
    << "DVyDy" << std::setw(spacing)
    << "swellingPressure" << std::setw(spacing)
    << "elasticEnergy" << std::setw(spacing)
    << "mixingEnergy" << std::setw(spacing)
    << "internalEnergy" << std::endl;

    for (int i=0; i<refArea.size();i++ ) {
        myfile
        <<  std::setprecision(15)
        << i << std::setw(spacing)
        << refArea[i] << std::setw(spacing)
        << areaRatio[i] << std::setw(spacing)
        << defGradXX[i] << std::setw(spacing)
        << defGradXY[i] << std::setw(spacing)
        << defGradYX[i] << std::setw(spacing)
        << defGradYY[i] << std::setw(spacing)
        << defGradXX[i]-prev_defGradXX[i] << std::setw(spacing)
        << defGradXY[i]-prev_defGradXY[i] << std::setw(spacing)
        << defGradYX[i]-prev_defGradYX[i] << std::setw(spacing)
        << defGradYY[i]-prev_defGradYY[i] << std::setw(spacing)
        << PK1stressXX[i] << std::setw(spacing)
        << PK1stressXY[i] << std::setw(spacing)
        << PK1stressYX[i] << std::setw(spacing)
        << PK1stressYY[i] << std::setw(spacing)
        << CstressXX[i] << std::setw(spacing)
        << CstressXY[i] << std::setw(spacing)
        << CstressYX[i] << std::setw(spacing)
        << CstressYY[i] << std::setw(spacing)
        << CstressXX[i]-prev_CstressXX[i] << std::setw(spacing)
        << CstressXY[i]-prev_CstressXY[i] << std::setw(spacing)
        << CstressYX[i]-prev_CstressYX[i] << std::setw(spacing)
        << CstressYY[i]-prev_CstressYY[i] << std::setw(spacing)
        << DtotVxDx[i] << std::setw(spacing)
        << DtotVxDy[i] << std::setw(spacing)
        << DtotVyDx[i] << std::setw(spacing)
        << DtotVyDy[i] << std::setw(spacing)
        << DhomVxDx[i] << std::setw(spacing)
        << DhomVxDy[i] << std::setw(spacing)
        << DhomVyDx[i] << std::setw(spacing)
        << DhomVyDy[i] << std::setw(spacing)
        << DVxDx[i] << std::setw(spacing)
        << DVxDy[i] << std::setw(spacing)
        << DVyDx[i] << std::setw(spacing)
        << DVyDy[i] << std::setw(spacing)
        << swellingPressurePerEle[i] << std::setw(spacing)
        << elasticEnergyPerEle[i] << std::setw(spacing)
        << mixingEnergyPerEle[i] << std::setw(spacing)
        << internalEnergyPerEle[i] << std::endl;
    }
    
    myfile.close();

}


void Configuration::dump_smoothcurves(const BaseSysData& baseData, const Parameters& pars, long& timeStep){
    
    int spacing =26;
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    if (pars.runMode == "stepShear" || pars.runMode == "surfaceShear"  || pars.runMode == "continuousShear"){
    
            if (pars.outputSubFolder=="none"){
                myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingStepNum))+"/smoothcurves-"+step+".txt").c_str());
            }else{
                myfile.open ((pars.outputFolderName +"/"+ pars.outputSubFolder+"/dsmoothcurves-"+step+".txt").c_str());
            }
    

    }else{
        myfile.open (pars.outputFolderName+"/smoothcurves-"+step+".txt");
    }
    myfile << "numMeshes:   " << baseData.numOriginalMeshes << std::endl;
    
    for (const auto& mesh: baseData.originalSurfaceMeshes)
    {
    
        for (const auto& node3: mesh.second)
        {
            int segment0 = baseData.nodeToSegments[node3][0];
            int segment1 = baseData.nodeToSegments[node3][1];
            
            int node2 =  baseData.surfaceSegments[segment0][0];
            double x2 = augmentedCurPosX[node2];
            double y2 = augmentedCurPosY[node2];
            
            int node4 =  baseData.surfaceSegments[segment1][1];
            double x4 = augmentedCurPosX[node4];
            double y4 = augmentedCurPosY[node4];
            
            double x3 = augmentedCurPosX[node3];
            double y3 = augmentedCurPosY[node3];
            
            
            
            for (double g =-1; g<=1; g+=0.05){
                
                struct point xy = HermitianInterpolation(x2, x3, x4, y2, y3, y4, pars.alpha_HermitPol, g);
                
                myfile
                <<  std::setprecision(15)
                << xy.x << std::setw(spacing)
                << xy.y << std::endl;
                
                
            }
        }
        myfile <<  "\n" << std::endl;
    }
    
    myfile.close();
    
    
    ////////////////
    std::ofstream myfile2;

    if (pars.runMode == "stepShear" || pars.runMode == "surfaceShear" || pars.runMode == "continuousShear"){
        
        if (pars.outputSubFolder=="none"){
                myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingStepNum))+"/gaps-"+step+".txt").c_str());
            }else{
                myfile.open ((pars.outputFolderName +"/"+ pars.outputSubFolder+"/gaps-"+step+".txt").c_str());
            }
    
    }else{
        myfile2.open (pars.outputFolderName+"/gaps-"+step+".txt");
    }

    
    for (const auto& g: interactions_nts)
    {
    
      
    myfile2
    <<  std::setprecision(15)
        << g.second[0] << std::setw(spacing)
        << g.second[1] << std::setw(spacing)
        << g.second[2] << std::setw(spacing)
        << g.second[3]<< std::endl;
    }
   myfile2 << "EOF"<< std::endl;
    
  
    
    myfile2.close();
    
    
}

