//
//#include <iostream>
//#include <string>
//#include <vector>
//#include <fstream>
//#include <sstream>
//#include <Eigen/Sparse>
//#include <Eigen/Dense>
//#include <cmath>
//#include <sys/stat.h>
//#include <dirent.h>
//#include <iomanip>
//#include <chrono>
//#include <valarray>
//#include "Parameters.hpp"
//#include "Configuration.hpp"
//#include "BaseSysData.hpp"
//
//
//void Configuration::compute_forces_walls(const BaseSysData& baseData, const Parameters& pars, const int& timeStep)
//{
//    auto start1 = std::chrono::high_resolution_clock::now();
//    contactsEnergy=0; //erase previous step data
//    //compute the deformation gradient
//    defGradXX = gradX * curPosX;
//    defGradYX = gradX * curPosY;
//    defGradXY = gradY * curPosX;
//    defGradYY = gradY * curPosY;
//    
////    auto finish1 = std::chrono::high_resolution_clock::now();
////    std::chrono::duration<double> elapsed1 = finish1 - start1;
////    std::cout << "elapsed time in defGrad:  " << elapsed1.count() << std::endl;
//    
//    //compute the determinant of defGrad
//    areaRatio = defGradXX.array() * defGradYY.array() - defGradXY.array() * defGradYX.array();
//    
//    
////    auto finish2 = std::chrono::high_resolution_clock::now();
////    std::chrono::duration<double> elapsed2 = finish2 - finish1;
////    std::cout << "elapsed time in areaRatio:  " << elapsed2.count() << std::endl;
//    
//    //compute the magnitude squared of defGrad
//    fSquared = defGradXX.array().pow(2)+defGradXY.array().pow(2)+defGradYX.array().pow(2)+defGradYY.array().pow(2);
//    
////    auto finish3 = std::chrono::high_resolution_clock::now();
////    std::chrono::duration<double> elapsed3 = finish3 - finish2;
////    std::cout << "elapsed time in fSquared:  " << elapsed3.count() << std::endl;
//    
//    //compute the local elastic energy density
//    elasticEnergyPerEle = pars.NkT/2.0 * (fSquared.array()- 2 - 2*log(abs(areaRatio.array())));
////    auto finish4 = std::chrono::high_resolution_clock::now();
////    std::chrono::duration<double> elapsed4 = finish4 - finish3;
////    std::cout << "elapsed time in elasticEnergyPerEle:  " << elapsed4.count() << std::endl;
////
//    //compute the local mixing energy density
//    mixingEnergyPerEle = pars.kTOverOmega*(areaRatio.array()-1.0)*(log((abs(areaRatio.array())-1.0)/abs(areaRatio.array()))+pars.chi/abs(areaRatio.array()));
////    auto finish5 = std::chrono::high_resolution_clock::now();
////    std::chrono::duration<double> elapsed5 = finish5 - finish4;
////    std::cout << "elapsed time in mixingEnergyPerEle:  " << elapsed5.count() << std::endl;
////
//    //compute the local total energy density
//    internalEnergyPerEle = elasticEnergyPerEle.array() + mixingEnergyPerEle.array();
////    auto finish6 = std::chrono::high_resolution_clock::now();
////    std::chrono::duration<double> elapsed6 = finish6 - finish5;
////    std::cout << "elapsed time in internalEnergyPerEle:  " << elapsed6.count() << std::endl;
//    
//    //compute the inverse of defGrad
//    invDefGradTransXX = defGradYY.array()/ areaRatio.array();
//    invDefGradTransXY = defGradXY.array()/ areaRatio.array() * -1 ;
//    invDefGradTransYX = defGradYX.array()/ areaRatio.array() * -1;
//    invDefGradTransYY = defGradXX.array()/ areaRatio.array();
//    
////    auto finish7 = std::chrono::high_resolution_clock::now();
////    std::chrono::duration<double> elapsed7 = finish7 - finish6;
////    std::cout << "elapsed time in invDefGradTrans:  " << elapsed7.count() << std::endl;
//    
//    swellingPressurePerEle = pars.kTOverOmega * ( (pars.chi+areaRatio.array()) / (areaRatio.array()).pow(2) + log((areaRatio.array()-1)/areaRatio.array()) );
//    
//    PK1stressXX = pars.NkT * defGradXX.array() + ( swellingPressurePerEle.array() * areaRatio.array() - pars.NkT ) * invDefGradTransXX.array();
//    PK1stressXY = pars.NkT * defGradXY.array() + ( swellingPressurePerEle.array() * areaRatio.array() - pars.NkT ) * invDefGradTransYX.array();
//    PK1stressYX = pars.NkT * defGradYX.array() + ( swellingPressurePerEle.array() * areaRatio.array() - pars.NkT ) * invDefGradTransXY.array();
//    PK1stressYY = pars.NkT * defGradYY.array() + ( swellingPressurePerEle.array() * areaRatio.array() - pars.NkT ) * invDefGradTransYY.array();
//    
//    CstressXX = pars.NkT/areaRatio.array()*(defGradXX.array()*defGradXX.array()+defGradXY.array()*defGradXY.array())+(swellingPressurePerEle.array()-1/areaRatio.array())*(defGradXX.array()*invDefGradTransXX.array()+defGradXY.array()*invDefGradTransYX.array());
//    CstressXY = pars.NkT/areaRatio.array()*(defGradXX.array()*defGradYX.array()+defGradXY.array()*defGradYY.array())+(swellingPressurePerEle.array()-1/areaRatio.array())*(defGradYX.array()*invDefGradTransXX.array()+defGradYY.array()*invDefGradTransYX.array());
//    CstressYX = pars.NkT/areaRatio.array()*(defGradYX.array()*defGradXX.array()+defGradYY.array()*defGradXY.array())+(swellingPressurePerEle.array()-1/areaRatio.array())*(defGradXX.array()*invDefGradTransXY.array()+defGradXY.array()*invDefGradTransYY.array());
//    CstressYY = pars.NkT/areaRatio.array()*(defGradYX.array()*defGradYX.array()+defGradYY.array()*defGradYY.array())+(swellingPressurePerEle.array()-1/areaRatio.array())*(defGradYX.array()*invDefGradTransXY.array()+defGradYY.array()*invDefGradTransYY.array());
//    
//    
//    //compute the nodal forces from the stresses
//    forceX = - gradX.transpose() * (PK1stressXX.array() * abs(refArea.array())).matrix() - gradY.transpose() * (PK1stressXY.array() * abs(refArea.array())).matrix();
//    forceY = - gradX.transpose() * (PK1stressYX.array() * abs(refArea.array())).matrix() - gradY.transpose() * (PK1stressYY.array() * abs(refArea.array())).matrix();
//    interForceX = forceX;
//    interForceY = forceY;
//    
//    
//    
//    if (pars.wallStyle=="harmonic"){
//        wallForceTop = pars.HWallStiffness * (0.5*(sign(curPosY.array()-topPos)+1))*(topPos-curPosY.array());
//        wallForceBottom = pars.HWallStiffness *(0.5*(sign(botPos-curPosY.array())+1))*(botPos-curPosY.array());
//        wallForceRight= pars.HWallStiffness * (0.5*(sign(curPosX.array()-rightPos)+1))*(rightPos-curPosX.array());
//        wallForceLeft= pars.HWallStiffness * (0.5*(sign(leftPos-curPosX.array())+1))*(leftPos-curPosX.array());
//        wallsEnergy = 0.5/pars.HWallStiffness * (wallForceTop.dot(wallForceTop) + wallForceBottom.dot(wallForceBottom) + wallForceLeft.dot(wallForceLeft) + wallForceRight.dot(wallForceRight));
//        
//        maxWallinterference = -(wallForceTop.minCoeff())/pars.HWallStiffness;
//        if(wallForceBottom.maxCoeff()/pars.HWallStiffness > maxWallinterference ){
//            maxWallinterference = wallForceBottom.maxCoeff()/pars.HWallStiffness;
//        }else if(-(wallForceRight.minCoeff())/pars.HWallStiffness > maxWallinterference ){
//            maxWallinterference = -(wallForceRight.minCoeff())/pars.HWallStiffness;
//        }else if (wallForceLeft.maxCoeff()/pars.HWallStiffness > maxWallinterference){
//            maxWallinterference = wallForceLeft.maxCoeff()/pars.HWallStiffness;
//        }
//    }
//    
//    if (pars.wallStyle=="powerlaw"){
//        wallForceTop=-(pars.PLWallEnergyScale/pars.PLWallLJScale)*12.0/((topPos-curPosY.array())/pars.PLWallLJScale).pow(13);
//        wallForceBottom=+(pars.PLWallEnergyScale/pars.PLWallLJScale)*12.0/((curPosY.array()-botPos)/pars.PLWallLJScale).pow(13);
//        
//        wallForceRight=-(pars.PLWallEnergyScale/pars.PLWallLJScale)*12.0/((rightPos-curPosX.array())/pars.PLWallLJScale).pow(13);
//        wallForceLeft=+(pars.PLWallEnergyScale/pars.PLWallLJScale)*12.0/((curPosX.array()-leftPos)/pars.PLWallLJScale).pow(13);
//        
//        wallsEnergy+=pars.PLWallEnergyScale*(1.0/(((topPos-curPosY.array())/pars.PLWallLJScale).pow(12)).sum());
//        wallsEnergy+=pars.PLWallEnergyScale*(1.0/(((curPosY.array()-botPos)/pars.PLWallLJScale).pow(12).sum()));
//        wallsEnergy+=pars.PLWallEnergyScale*(1.0/(((rightPos-curPosX.array())/pars.PLWallLJScale).pow(12)).sum());
//        wallsEnergy+=pars.PLWallEnergyScale*(1.0/(((curPosX.array()-rightPos)/pars.PLWallLJScale).pow(12)).sum());
//    }
//    
//    
//    forceX += wallForceLeft + wallForceRight;
//    forceY += wallForceBottom + wallForceTop;
//    
//    auto finish8 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed8 = finish8 - start1;
//    std::cout << "elapsed time in internal nodes:  " << elapsed8.count() << std::endl;
//    
//    
//    maxInterference = 0;
//    segmentIinteractions = 0;
//    nodeIinteractions = 0;
//    gaps.clear();
//    
//    numXBins = floor(lxNew*(1+2*pars.imagesMargin)/pars.verletCellCutoff);
//    numYBins = floor(lyNew*(1+2*pars.imagesMargin)/pars.verletCellCutoff);
//    verletCellSizeX  = lxNew/numXBins;
//    verletCellSizeY = lyNew/numYBins;
//    cellList.resize(numXBins*numYBins+baseData.numSurfaceNodes,-1);
//    
//    update_cells(baseData, pars);
//    
//    for (int m1y=0; m1y<numYBins; m1y++)
//    {
//        
//        for (int m1x=0; m1x<numXBins; m1x++)
//        {
//            int m1 = numXBins*m1y + m1x + baseData.numSurfaceNodes;
//            for (auto const& delta: neighborBinDelta)
//            {
//                int  nBinXid = m1x + delta.first;
//                int  nBinYid = m1y + delta.second;
//                
//                if ((nBinXid < 0) || (nBinYid < 0) || (nBinXid == numXBins) || (nBinYid == numYBins))
//                {
//                    continue;
//                    
//                }
//                int   m2 = numXBins*nBinYid + nBinXid + baseData.numSurfaceNodes;
//                int slaveNodeId = cellList[m1];  //this is the slave surface node id in flatSurfaceNodes vector, not the global id
//                while (slaveNodeId>=0)
//                {
//                    int slaveMesh = baseData.nodeToSegments[baseData.flatSurfaceNodes[slaveNodeId]][2];
//                    int masterNodeId = cellList[m2]; //this is the master surface node id in flatSurfaceNodes vector, not the global id
//                    while (masterNodeId>=0)
//                    {
//                        int masterMesh = baseData.nodeToSegments[baseData.flatSurfaceNodes[masterNodeId]][2];
//                        if (slaveMesh==masterMesh)
//                        {
//                            masterNodeId = cellList[masterNodeId];
//                            continue;
//                        }
////                        std::cout << " sNode  " << baseData.flatSurfaceNodes[slaveNodeId] << std::endl;
////                        std::cout << " mNode  " << baseData.flatSurfaceNodes[masterNodeId] << std::endl<< std::endl;
//                        int segment0 = baseData.nodeToSegments[baseData.flatSurfaceNodes[masterNodeId]][0];
//                        int segment1 = baseData.nodeToSegments[baseData.flatSurfaceNodes[masterNodeId]][1];
//                        NTS_interaction(baseData.flatSurfaceNodes[slaveNodeId],segment0, baseData, pars);
//                        NTS_interaction(baseData.flatSurfaceNodes[slaveNodeId],segment1, baseData, pars);
//                        
//                        masterNodeId = cellList[masterNodeId];
//                    }
//                    
//                    slaveNodeId = cellList[slaveNodeId];
//                }
//            }
//        }
//    }
//    
//    auto finish9 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed9 = finish9 - finish8;
//    std::cout << "elapsed time in new verlet cells:  " << elapsed9.count() << std::endl;
//    
//    
//    for (auto const& nodeRow: gaps)
//    {
//        
//        std::vector<double> shortestPath = nodeRow.second;
//        
//        if (shortestPath[1]>= 0)
//        {
//            continue;
//        }
//        if (shortestPath[0] > maxInterference)
//        {
//            maxInterference = shortestPath[0];
//        }
//        if (shortestPath.size()==20)
//        {
//            
//            segmentIinteractions++ ;
//            if (nodeRow.first.first < baseData.numOriginalNodes)
//            {
//                forceX(nodeRow.first.first) = forceX(nodeRow.first.first) + shortestPath[4] ;
//                forceY(nodeRow.first.first) = forceY(nodeRow.first.first) + shortestPath[5] ;
//            }
//            if ( shortestPath[6] < baseData.numOriginalNodes){
//                forceX(shortestPath[6]) = forceX(shortestPath[6])  + shortestPath[8] ;
//                forceY(shortestPath[6]) = forceY(shortestPath[6]) + shortestPath[9];
//            }
//            if ( shortestPath[10] < baseData.numOriginalNodes){
//                forceX(shortestPath[10]) = forceX(shortestPath[10]) + shortestPath[12];
//                forceY(shortestPath[10]) = forceY(shortestPath[10]) + shortestPath[13];
//            }
//            
//            
//            contactsEnergy += pars.penaltyStiffness/2 *(shortestPath[0]*shortestPath[0]);
//            
//            
//        }else{
//            
//            nodeIinteractions++;
//            if (nodeRow.first.first < baseData.numOriginalNodes){
//                forceX(nodeRow.first.first) = forceX(nodeRow.first.first) + shortestPath[4] ;
//                forceY(nodeRow.first.first) = forceY(nodeRow.first.first) + shortestPath[5] ;
//            }
//            if ( shortestPath[6] < baseData.numOriginalNodes){
//                forceX(shortestPath[6]) = forceX(shortestPath[6])  + shortestPath[8] ;
//                forceY(shortestPath[6]) = forceY(shortestPath[6]) + shortestPath[9];
//            }
//            
//            contactsEnergy += pars.penaltyStiffness/2 *(shortestPath[0]*shortestPath[0]);
//            
//        }
//    }
//    
//    
//    internalEnergy = internalEnergyPerEle.dot(refArea) + contactsEnergy;
//    totalEnergy= internalEnergy + wallsEnergy;
//    auto finish10 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed10 = finish10 - finish9;
//    std::cout << "elapsed time in gaps map looping:  " << elapsed10.count() << std::endl;
//    
//    std::cout << "segmentIinteractions  " << segmentIinteractions << std::endl;
//    std::cout << "nodeIinteractions  " << nodeIinteractions << std::endl;
//    std::cout << "max skin interference   " << maxInterference << std::endl;
//    std::cout << "max wall interference   " << maxWallinterference << std::endl;
//    std::cout << "min J  " <<areaRatio.array().minCoeff() << std::endl;
//    
//}
