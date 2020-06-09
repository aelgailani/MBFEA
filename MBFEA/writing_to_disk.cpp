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

void Configuration::dump_facets(const BaseSysData& baseData, const Parameters& pars, long& timeStep){
    
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    if (pars.runMode == "stepShear" || pars.runMode == "continuousShear"){
    myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingTimeStep))+"/facets-"+step+".txt").c_str());
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


void Configuration::dump_global_data(const Parameters& pars, const long& timeStep, std::string mode, std::string purpose){
    int spacing =26;
    std::string fname ;
    std::string first = std::to_string(timeStep/pars.splitDataEvery*pars.splitDataEvery);
    std::string last =std::to_string(timeStep/pars.splitDataEvery*pars.splitDataEvery+pars.splitDataEvery-1);
    
    
    if (purpose=="running"){
        fname = "/data-steps-"+first+"-"+last+".txt";
    }else if (purpose=="final"){
        fname = "/final-data.txt";
    }
   
    if (mode=="write" || first != lastStepFirst){
        std::ofstream dataFile;
        if (pars.runMode == "stepShear" || pars.runMode == "continuousShear"){
        dataFile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingTimeStep))).c_str()+fname);
        }else{
            dataFile.open (pars.outputFolderName+fname);
        }
        dataFile
        <<  "step"  << std::setw(25)
        <<  "totalEnergy"  << std::setw(spacing)
        <<  "contactEnergy"  << std::setw(spacing)
        <<  "shearVirial"  << std::setw(spacing)
        <<  "pressureVirial"  << std::setw(spacing)
        <<  "maxResidualF"  << std::setw(spacing)
        <<  "pressure1"  << std::setw(spacing)
        <<  "pressure2"  << std::setw(spacing)
        <<  "KWpressure2"  << std::setw(spacing)
        <<  "shearStress1"  << std::setw(spacing)
        <<  "shearStress2"  << std::setw(spacing)
        <<  "KWshearStress2"  << std::setw(spacing)
        <<  "boxArea"  << std::setw(spacing)
        <<  "materialArea"  << std::setw(spacing)
        <<  "phi"  << std::setw(spacing)
        <<  "e0"  << std::setw(spacing)
        <<  "e1" << std::setw(spacing)
        <<  "DPOverDe0" << std::setw(spacing)
        <<  "DSOverDe1qq" << std::setw(spacing)
        <<  "dt" << std::setw(spacing)
        <<  "defRate" << std::setw(spacing)
        <<  "penalty" << std::endl;
        dataFile.close();
        
        
    }else if (mode=="append"){
        std::ofstream dataFile;
        if (pars.runMode == "stepShear" || pars.runMode == "continuousShear"){
        dataFile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingTimeStep))).c_str()+fname,  std::ios_base::app);
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
        <<  P1  << std::setw(spacing)
        <<  P2  << std::setw(spacing)
        <<  (KWoodYY+KWoodXX)/(2*LX*LY)  << std::setw(spacing)
        <<  S1  << std::setw(spacing)
        <<  S2  << std::setw(spacing)
        <<  (KWoodYY-KWoodXX)/(2*LX*LY)   << std::setw(spacing)
        <<  A  << std::setw(spacing)
        <<  A_material  << std::setw(spacing)
        <<  phi  << std::setw(spacing)
        <<  e0  << std::setw(spacing)
        <<  e1 << std::setw(spacing)
        <<  DPOverDe0 << std::setw(spacing)
        <<  DSOverDe1 << std::setw(spacing)
        <<  pars.dt  << std::setw(spacing)
        <<  pars.deformationRate  << std::setw(spacing)
        <<  pars.ntsHarmonicPenaltyStiffness << std::endl;
        dataFile.close();
        
        
    }
    
   lastStepFirst = first;
}


void Configuration::dump_per_node(const BaseSysData& baseData, const Parameters& pars, long& timeStep){
    int spacing =26;
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    if (pars.runMode == "stepShear" || pars.runMode == "continuousShear"){
    myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingTimeStep))+"/data-per-node-"+step+".txt").c_str());
    }else{
        myfile.open (pars.outputFolderName+"/data-per-node-"+step+".txt");
    }
    
    myfile << "Basic_data:" << std::endl;
    myfile << "timeStep" << "\t" << timeStep << std::endl;
    myfile << "number_of_nodes" << "\t" << baseData.numOriginalNodes << std::endl;
    myfile << "kTOverOmega" << "\t" << pars.kTOverOmega << std::endl;
    myfile << "phi" << "\t" << phi << std::endl;
    myfile << "e0_strain" << "\t" << e0 << std::endl;
    myfile << "e1_strain" << "\t" << e1 << std::endl;
    myfile << "dt" << "\t" << pars.dt << std::endl;
    myfile << "max_residual" << "\t" << maxR << std::endl;
    myfile << "deformation_rate" << "\t" << pars.deformationRate   << std::endl;
    myfile << "penalty_stiffness" << "\t" << pars.ntsHarmonicPenaltyStiffness << std::endl;
    myfile << "verlet_cell_cutoff" << "\t" << pars.verletCellCutoff << std::endl;
    myfile << "original_box_LRBT" << "\t" << leftPos << "\t" << rightPos << "\t" << botPos << "\t" << topPos<< "\n"  << std::endl;
    myfile << "Nodes_data:" << std::endl;
    myfile
    << "id"  << std::setw(spacing)
    <<  "x"  << std::setw(spacing)
    <<  "y"  << std::setw(spacing)
    <<  "fx"  << std::setw(spacing)
    <<  "fy"  << std::setw(spacing)
    <<  "contactFx"  << std::setw(spacing)
    <<  "contactFy"  << std::endl;
    
    for (int i=0; i<curPosX.size();i++ ) {
        myfile
        <<  std::setprecision(15)
        << i << std::setw(spacing)
        << curPosX[i] << std::setw(spacing)
        << curPosY[i] << std::setw(spacing)
        << forceX[i] << std::setw(spacing)
        << forceY[i] << std::setw(spacing)
        << surfaceForceX[i] << std::setw(spacing)
        << surfaceForceY[i] << std::endl;
    }
    
    myfile.close();

}

void Configuration::dump_per_node_periodic_images_on(const BaseSysData& baseData, const Parameters& pars, long& timeStep){
    int spacing =26;
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    if (pars.runMode == "stepShear" || pars.runMode == "continuousShear"){
    myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingTimeStep))+"/data-per-node-periodic-"+step+".txt").c_str());
    }else{
        myfile.open (pars.outputFolderName+"/data-per-node-periodic-"+step+".txt");
    }
    myfile << "Basic_data:" << std::endl;
    myfile << "timeStep" << "\t" << timeStep << std::endl;
    myfile << "number_of_augmented_nodes" << "\t" << baseData.numNodes << std::endl;
    myfile << "images_margin" << "\t" << pars.imagesMargin << std::endl;
    myfile << "effective_simulation_box_left_boundary" << "\t" << leftPos - pars.imagesMargin*lxNew << std::endl;
    myfile << "effective_simulation_box_bot_boundary" << "\t" << botPos - pars.imagesMargin*lyNew << std::endl;
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
    DVxDx = gradX * forceX;
    DVxDy = gradY * forceX;
    DVyDx = gradY * forceY;
    DVyDy = gradY * forceY;
    
    
    std::string step = std::to_string(timeStep);
    std::ofstream myfile;
    if (pars.runMode == "stepShear" || pars.runMode == "continuousShear"){
    myfile.open ((pars.outputFolderName +"/step-"+std::to_string(int(pars.startingTimeStep))+"/data-per-ele-"+step+".txt").c_str());
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
    << "PK1StressXX" << std::setw(spacing)
    << "PK1StressXY" << std::setw(spacing)
    << "PK1StressYX" << std::setw(spacing)
    << "PK1StressYY" << std::setw(spacing)
    << "CStressXX" << std::setw(spacing)
    << "CStressXY" << std::setw(spacing)
    << "CStressYX" << std::setw(spacing)
    << "CStressYY" << std::setw(spacing)
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
        << PK1stressXX[i] << std::setw(spacing)
        << PK1stressXY[i] << std::setw(spacing)
        << PK1stressYX[i] << std::setw(spacing)
        << PK1stressYY[i] << std::setw(spacing)
        << CstressXX[i] << std::setw(spacing)
        << CstressXY[i] << std::setw(spacing)
        << CstressYX[i] << std::setw(spacing)
        << CstressYY[i] << std::setw(spacing)
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
