
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
#include "Parameters.hpp"
#include "Configuration.hpp"
#include "BaseSysData.hpp"


void Configuration::compute_forces_walls(const BaseSysData& baseData, const Parameters& pars)
{
    
    contactsEnergy=0; //erase previous step data

    //compute the deformation gradient
    defGradXX = gradX * curPosX;
    defGradYX = gradX * curPosY;
    defGradXY = gradY * curPosX;
    defGradYY = gradY * curPosY;
    
    //compute the determinant of defGrad
    areaRatio = defGradXX.array() * defGradYY.array() - defGradXY.array() * defGradYX.array();
    std::cout << "min J  " <<areaRatio.array().minCoeff() << std::endl;
    
    //compute the magnitude squared of defGrad
    fSquared = defGradXX.array().pow(2)+defGradXY.array().pow(2)+defGradYX.array().pow(2)+defGradYY.array().pow(2);
    
    //compute the local elastic energy density
    elasticEnergyPerEle = pars.NkT/2.0 * (fSquared.array()- 2 - 2*log(abs(areaRatio.array())));
    
    //compute the local mixing energy density
    mixingEnergyPerEle = pars.kTOverOmega*(areaRatio.array()-1.0)*(log((abs(areaRatio.array())-1.0)/abs(areaRatio.array()))+pars.chi/abs(areaRatio.array()));
    
    //compute the local total energy density
    internalEnergyPerEle = elasticEnergyPerEle.array() + mixingEnergyPerEle.array();
    
    //compute the inverse of defGrad
    invDefGradTransXX = defGradYY.array()/ areaRatio.array();
    invDefGradTransXY = defGradXY.array()/ areaRatio.array() * -1 ;
    invDefGradTransYX = defGradYX.array()/ areaRatio.array() * -1;
    invDefGradTransYY = defGradXX.array()/ areaRatio.array();
    
    swellingPressurePerEle = pars.kTOverOmega * ( (pars.chi+areaRatio.array()) / (areaRatio.array()).pow(2) + log((areaRatio.array()-1)/areaRatio.array()) );
    
    PK1stressXX = pars.NkT * defGradXX.array() + ( swellingPressurePerEle.array() * areaRatio.array() - pars.NkT ) * invDefGradTransXX.array();
    PK1stressXY = pars.NkT * defGradXY.array() + ( swellingPressurePerEle.array() * areaRatio.array() - pars.NkT ) * invDefGradTransYX.array();
    PK1stressYX = pars.NkT * defGradYX.array() + ( swellingPressurePerEle.array() * areaRatio.array() - pars.NkT ) * invDefGradTransXY.array();
    PK1stressYY = pars.NkT * defGradYY.array() + ( swellingPressurePerEle.array() * areaRatio.array() - pars.NkT ) * invDefGradTransYY.array();
    
    CstressXX = pars.NkT/areaRatio.array()*(defGradXX.array()*defGradXX.array()+defGradXY.array()*defGradXY.array())+(swellingPressurePerEle.array()-1/areaRatio.array())*(defGradXX.array()*invDefGradTransXX.array()+defGradXY.array()*invDefGradTransYX.array());
    CstressXY = pars.NkT/areaRatio.array()*(defGradXX.array()*defGradYX.array()+defGradXY.array()*defGradYY.array())+(swellingPressurePerEle.array()-1/areaRatio.array())*(defGradYX.array()*invDefGradTransXX.array()+defGradYY.array()*invDefGradTransYX.array());
    CstressYX = pars.NkT/areaRatio.array()*(defGradYX.array()*defGradXX.array()+defGradYY.array()*defGradXY.array())+(swellingPressurePerEle.array()-1/areaRatio.array())*(defGradXX.array()*invDefGradTransXY.array()+defGradXY.array()*invDefGradTransYY.array());
    CstressYY = pars.NkT/areaRatio.array()*(defGradYX.array()*defGradYX.array()+defGradYY.array()*defGradYY.array())+(swellingPressurePerEle.array()-1/areaRatio.array())*(defGradYX.array()*invDefGradTransXY.array()+defGradYY.array()*invDefGradTransYY.array());
    
    
    
    
    
    //compute the nodal forces from the stresses
    forceX = - gradX.transpose() * (PK1stressXX.array() * abs(refArea.array())).matrix() - gradY.transpose() * (PK1stressXY.array() * abs(refArea.array())).matrix();
    forceY = - gradX.transpose() * (PK1stressYX.array() * abs(refArea.array())).matrix() - gradY.transpose() * (PK1stressYY.array() * abs(refArea.array())).matrix();
    interForceX = forceX;
    interForceY = forceY;
    
    double maxWallinterference = 0;
    //do harmonic repulsion on any nodes above topPos or below botPos
    if (pars.wallStyle=="harmonic"){
        wallForceTop = pars.HWallStiffness * (0.5*(sign(curPosY.array()-topPos)+1))*(topPos-curPosY.array());
        wallForceBottom = pars.HWallStiffness *(0.5*(sign(botPos-curPosY.array())+1))*(botPos-curPosY.array());
        wallForceRight= pars.HWallStiffness * (0.5*(sign(curPosX.array()-rightPos)+1))*(rightPos-curPosX.array());
        wallForceLeft= pars.HWallStiffness * (0.5*(sign(leftPos-curPosX.array())+1))*(leftPos-curPosX.array());
        wallsEnergy = 0.5/pars.HWallStiffness * (wallForceTop.dot(wallForceTop) + wallForceBottom.dot(wallForceBottom) + wallForceLeft.dot(wallForceLeft) + wallForceRight.dot(wallForceRight));

        maxWallinterference = -(wallForceTop.minCoeff())/pars.HWallStiffness;
        if(wallForceBottom.maxCoeff()/pars.HWallStiffness > maxWallinterference ){
            maxWallinterference = wallForceBottom.maxCoeff()/pars.HWallStiffness;
        }else if(-(wallForceRight.minCoeff())/pars.HWallStiffness > maxWallinterference ){
            maxWallinterference = -(wallForceRight.minCoeff())/pars.HWallStiffness;
        }else if (wallForceLeft.maxCoeff()/pars.HWallStiffness > maxWallinterference){
            maxWallinterference = wallForceLeft.maxCoeff()/pars.HWallStiffness;
        }
    }
    
    if (pars.wallStyle=="powerlaw"){
        wallForceTop=-(pars.PLWallEnergyScale/pars.PLWallLJScale)*12.0/((topPos-curPosY.array())/pars.PLWallLJScale).pow(13);
        wallForceBottom=+(pars.PLWallEnergyScale/pars.PLWallLJScale)*12.0/((curPosY.array()-botPos)/pars.PLWallLJScale).pow(13);
        
        wallForceRight=-(pars.PLWallEnergyScale/pars.PLWallLJScale)*12.0/((rightPos-curPosX.array())/pars.PLWallLJScale).pow(13);
        wallForceLeft=+(pars.PLWallEnergyScale/pars.PLWallLJScale)*12.0/((curPosX.array()-leftPos)/pars.PLWallLJScale).pow(13);
        
        wallsEnergy+=pars.PLWallEnergyScale*(1.0/(((topPos-curPosY.array())/pars.PLWallLJScale).pow(12)).sum());
        wallsEnergy+=pars.PLWallEnergyScale*(1.0/(((curPosY.array()-botPos)/pars.PLWallLJScale).pow(12).sum()));
        wallsEnergy+=pars.PLWallEnergyScale*(1.0/(((rightPos-curPosX.array())/pars.PLWallLJScale).pow(12)).sum());
        wallsEnergy+=pars.PLWallEnergyScale*(1.0/(((curPosX.array()-rightPos)/pars.PLWallLJScale).pow(12)).sum());
    }
    
    
    forceX += wallForceLeft + wallForceRight;
    forceY += wallForceBottom + wallForceTop;
    
    
    //compute surface forces
    ////construct and fill the bins
    double LX = ((curPosX.maxCoeff())+0.1)-((curPosX.minCoeff())-0.1);
    double LY = ((curPosY.maxCoeff())+0.1)-((curPosY.minCoeff())-0.1);
    double numXBins = floor(LX/pars.verletCellCutoff);
    double numYBins = floor(LY/pars.verletCellCutoff);
    double verletCellSizeX = LX/numXBins;
    double verletCellSizeY = LY/numYBins;
    double farLeft = (curPosX.minCoeff()) - 0.1;
    double farBottom = (curPosY.minCoeff()) - 0.1;
    int xBin, yBin;
    std::map<std::pair<int,int>, std::vector<int>> spatialGridNodes,spatialGridSegments,spatialGridMeshes;
    std::map<std::pair<int,int>, std::vector<double>> gaps;
    
    
    for (int meshID=0; meshID < baseData.numMeshes; meshID++)
    {
        for (int nodeID=0; nodeID < baseData.numNodesPerMesh; nodeID++)
        {
            xBin = int( floor( (curPosX[baseData.surfaceMeshes.at(meshID)[nodeID]] - farLeft)/verletCellSizeX) );
            yBin = int( floor( (curPosY[baseData.surfaceMeshes.at(meshID)[nodeID]] - farBottom)/verletCellSizeY) );
            spatialGridNodes[{xBin,yBin}].push_back(baseData.surfaceMeshes.at(meshID)[nodeID]);
            spatialGridSegments[{xBin,yBin}].push_back(baseData.nodeToSegments[baseData.surfaceMeshes.at(meshID)[nodeID]][0]);  //add the first segment of this particle
            spatialGridSegments[{xBin,yBin}].push_back(baseData.nodeToSegments[baseData.surfaceMeshes.at(meshID)[nodeID]][1]);  //add the second segment of this particle
            spatialGridMeshes[{xBin,yBin}].push_back(meshID); //add the mesh of this particle
            
        }
    }
    
    double maxInterference = 0;
    
    
    //// loop over the bins to compute forces
    std::pair<int,int> neighborBinDelta[9] = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,0},{0,1},{1,-1},{1,0},{1,1}};
    
    for (auto const& bin : spatialGridNodes)
    {
        if (bin.second.size()==0) // check if the bin is empty
        {
            continue;
        }
        
        std::vector<int> segments = spatialGridSegments[bin.first];
        std::vector<int> meshes = spatialGridMeshes[bin.first];
        unsigned long meshesNum = meshes.size();
        
        for (auto const& delta: neighborBinDelta)
        {
            std::pair<int,int> neighborBinKey = {bin.first.first + delta.first, bin.first.second + delta.second};
            
            if (spatialGridSegments.find(neighborBinKey) == spatialGridSegments.end()) // check if neighbour bin key is nonexistant
            {
                continue;
            }
            
            if (spatialGridMeshes[neighborBinKey].size()<1) // check if the neighbour bin is empty
            {
                continue;
            }
            
            for (int const& i : spatialGridSegments[neighborBinKey])
            {
                if(std::find(segments.begin(), segments.end(), i) != segments.end()) {
                    continue;
                } else {
                    segments.push_back(i);
                }
            }
            
            for (int const& i : spatialGridMeshes[neighborBinKey])
            {
                if(std::find(meshes.begin(), meshes.end(), i) != meshes.end()) {
                    continue;
                } else {
                    meshes.push_back(i);
                }
            }
            meshesNum = meshes.size();
            
        }
        
        if  (baseData.numMeshes<2) // check if all nodes live in the same mesh
        {
            continue;
        }
        
        for (int const& node: bin.second)
        {
            int slaveMesh = baseData.nodeToSegments[node][2];
            for(int const& segment: segments)
            {
                int masterMesh = baseData.surfaceSegments[segment][2];
                
                if (slaveMesh==masterMesh)
                {
                    continue;
                }
                
                int node0 = baseData.surfaceSegments[segment][0];
                int node1 = baseData.surfaceSegments[segment][1];
                
                double xi = curPosX[node];
                double yi = curPosY[node];
                double x0 = curPosX[node0];
                double y0 = curPosY[node0];
                double x1 = curPosX[node1];
                double y1 = curPosY[node1];
                
                //exclude segments that lies outside the walls to solve the sliiping failure at high phi
                if ((y0>topPos && y1>topPos) || (y0<botPos && y1<botPos) || (x0<leftPos && x1<leftPos)  || (x0>rightPos && x1>rightPos)){
                    continue;
                }
                
                
                
                double dx = x1-x0;
                double dy = y1-y0;
                double L = sqrt(std::pow(dx,2)+std::pow(dy,2));
                double nx = dy/L;
                double ny = - dx/L;
                double s = (xi-x0)*dx/std::pow(L,2)+(yi-y0)*dy/std::pow(L,2);
                double gap = (xi-x0)*nx+(yi-y0)*ny;
                double gapSign = 1;
                if (gap<0){ gapSign = -1; }
                
                if (s>=0 && s<=1)
                {
                    double f = pars.penaltyStiffness * abs(gap);
                    double f0 = -f*(1-s);
                    double f1 = -f * s;
                    double fx = f * (nx);
                    double fy = f * (ny);
                    double f0x = -fx * (1-s);
                    double f0y = -fy * (1-s);
                    double f1x = -fx * (s);
                    double f1y = -fy * (s);
                    
                    if (gaps[{node,masterMesh}].size()==0) {
                        gaps[{node,masterMesh}]={abs(gap),gapSign,double(node),f,fx,fy,double(node0),f0,f0x,f0y,double(node1),f1,f1x,f1y,nx,ny,s,double(segment),(xi-x0)*nx,(yi-y0)*ny };
                        continue;
                    }
                    if (gaps[{node,masterMesh}][0]>abs(gap)){
                        gaps[{node,masterMesh}] = {abs(gap),gapSign,double(node),f,fx,fy,double(node0),f0,f0x,f0y,double(node1),f1,f1x,f1y,nx,ny,s,double(segment),(xi-x0)*nx,(yi-y0)*ny };
                        continue;
                    };
                }else{
                    
                    double r0ix = xi-x0 ;
                    double r0iy = yi-y0;
                    double r1ix = xi-x1;
                    double r1iy = yi-y1;
                    
                    double gap0 = sqrt( std::pow(r0ix,2) + std::pow(r0iy,2) );
                    double gap1 = sqrt( std::pow(r1ix,2) + std::pow(r1iy,2) );
                    
                    double f0i= pars.penaltyStiffness  * gap0;
                    double f0ix = -pars.penaltyStiffness  * r0ix;
                    double f0iy = -pars.penaltyStiffness  * r0iy;
                    double f00= -f0i;
                    double f00x = -f0ix;
                    double f00y = -f0iy;
                    
                    double f1i = pars.penaltyStiffness  * gap1;
                    double f1ix = -pars.penaltyStiffness  * r1ix;
                    double f1iy = -pars.penaltyStiffness  * r1iy;
                    double f11 = -f1i;
                    double f11x = -f1ix;
                    double f11y = -f1iy;
                    
                    if (gaps[{node,masterMesh}].size()==0){
                        gaps[{node,masterMesh}] = {gap0,gapSign,double(node),f0i,f0ix,f0iy,double(node0),f00,f00x,f00y,double(segment),r0ix,r0iy};
                        
                    }
                    
                    if (gaps[{node,masterMesh}][0]>gap0){
                        gaps[{node,masterMesh}] = {gap0,gapSign,double(node),f0i,f0ix,f0iy,double(node0),f00,f00x,f00y,double(segment),r0ix,r0iy};
                        
                    }
                    
                    if (gaps[{node,masterMesh}][0]>gap1){
                        gaps[{node,masterMesh}] = {gap1,gapSign,double(node),f1i,f1ix,f1iy,double(node1),f11,f11x,f11y,double(segment),r1ix,r1iy};
                        
                    }
                }
                
            }
        }
        
    }
    
    int segmentIinteractions = 0;
    int nodeIinteractions = 0;
    
    for (auto const& nodeRow: gaps)
    {
        
        std::vector<double> shortestPath = nodeRow.second;
        
        if (shortestPath[1]>= 0) {continue;}
        if (shortestPath.size()==20)
        {
            if (shortestPath[0] > maxInterference){
                maxInterference = shortestPath[0];
            }
            segmentIinteractions++ ;
            if (nodeRow.first.first < baseData.numOriginalNodes){
                forceX(nodeRow.first.first) = forceX(nodeRow.first.first) + shortestPath[4] ;
                forceY(nodeRow.first.first) = forceY(nodeRow.first.first) + shortestPath[5] ;
            }
            if ( shortestPath[6] < baseData.numOriginalNodes){
                forceX(shortestPath[6]) = forceX(shortestPath[6])  + shortestPath[8] ;
                forceY(shortestPath[6]) = forceY(shortestPath[6]) + shortestPath[9];
            }
            if ( shortestPath[10] < baseData.numOriginalNodes){
                forceX(shortestPath[10]) = forceX(shortestPath[10]) + shortestPath[12];
                forceY(shortestPath[10]) = forceY(shortestPath[10]) + shortestPath[13];
            }
         
            
            contactsEnergy += pars.penaltyStiffness/2 *(shortestPath[0]*shortestPath[0]);

        }else{
            
            nodeIinteractions++;
            if (nodeRow.first.first < baseData.numOriginalNodes){
                forceX(nodeRow.first.first) = forceX(nodeRow.first.first) + shortestPath[4] ;
                forceY(nodeRow.first.first) = forceY(nodeRow.first.first) + shortestPath[5] ;
            }
            if ( shortestPath[6] < baseData.numOriginalNodes){
                forceX(shortestPath[6]) = forceX(shortestPath[6])  + shortestPath[8] ;
                forceY(shortestPath[6]) = forceY(shortestPath[6]) + shortestPath[9];
            }
            
            contactsEnergy += pars.penaltyStiffness/2 *(shortestPath[0]*shortestPath[0]);

        }
        
    }
    
    internalEnergy = internalEnergyPerEle.dot(refArea) + contactsEnergy;
    totalEnergy= internalEnergy + wallsEnergy;
    
    std::cout << "segmentIinteractions  " << segmentIinteractions << std::endl;
    std::cout << "nodeIinteractions  " << nodeIinteractions << std::endl;
    std::cout << "max skin interference   " << maxInterference << std::endl;
    std::cout << "max wall interference   " << maxWallinterference << std::endl;

}
