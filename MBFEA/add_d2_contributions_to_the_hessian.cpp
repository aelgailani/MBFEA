
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
#include <stdio.h>

void Configuration::add_d2_contributions_to_the_hessian(double penaltyStifness, double xs,double ys,double xm,double ym, int snode, int mnode, const BaseSysData& baseData){
    
    double d2 = pow((pow((xs - xm),2) + pow((ys - ym),2)),0.5);
    double Dd2Dxs = (1.*(-xm + xs))/pow(pow(xm - xs,2) + pow(ym - ys,2),0.5);
    double Dd2Dys = (1.*(-ym + ys))/pow(pow(xm - xs,2) + pow(ym - ys,2),0.5);
    double Dd2Dxm = (-1.*(-xm + xs))/pow(pow(xm - xs,2) + pow(ym - ys,2),0.5);
    double Dd2Dym = (-1.*(-ym + ys))/pow(pow(xm - xs,2) + pow(ym - ys,2),0.5);
    
    double DDd2DxsDxs = (-1.*pow(xm - 1.*xs,2))/pow(pow(xm - xs,2) + pow(ym - ys,2),1.5) + 1./pow(pow(xm - xs,2) + pow(ym - ys,2),0.5);
    double DDd2DxsDys = (-1.*(xm - 1.*xs)*(ym - 1.*ys))/pow(pow(xm - xs,2) + pow(ym - ys,2),1.5);
    double DDd2DxsDxm = (-1.*(xm - xs)*(-xm + xs))/pow(pow(xm - xs,2) + pow(ym - ys,2),1.5) - 1./pow(pow(xm - xs,2) + pow(ym - ys,2),0.5);
    double DDd2DxsDym = (-1.*(-xm + xs)*(ym - ys))/pow(pow(xm - xs,2) + pow(ym - ys,2),1.5);
    
    double DDd2DysDys = 1./pow(pow(xm - xs,2) + pow(ym - ys,2),0.5) - (1.*pow(ym - 1.*ys,2))/pow(pow(xm - xs,2) + pow(ym - ys,2),1.5);
    double DDd2DysDxm = (1.*(xm - 1.*xs)*(ym - 1.*ys))/pow(pow(xm - xs,2) + pow(ym - ys,2),1.5);
    double DDd2DysDym = -1./pow(pow(xm - xs,2) + pow(ym - ys,2),0.5) + (1.*pow(ym - 1.*ys,2))/pow(pow(xm - xs,2) + pow(ym - ys,2),1.5);
    
    double DDd2DxmDxm = (-1.*pow(xm - 1.*xs,2))/pow(pow(xm - xs,2) + pow(ym - ys,2),1.5) + 1./pow(pow(xm - xs,2) + pow(ym - ys,2),0.5);
    double DDd2DxmDym = (-1.*(xm - 1.*xs)*(ym - 1.*ys))/pow(pow(xm - xs,2) + pow(ym - ys,2),1.5);
    
    double DDd2DymDym = 1./pow(pow(xm - xs,2) + pow(ym - ys,2),0.5) - (1.*pow(ym - 1.*ys,2))/pow(pow(xm - xs,2) + pow(ym - ys,2),1.5);
    
    
    
    
     if ( snode < baseData.numOriginalNodes && mnode < baseData.numOriginalNodes){
           
         Hixjx.coeffRef(snode,mnode)+= penaltyStifness *(Dd2Dxs*Dd2Dxm+d2*DDd2DxsDxm);
         Hixjy.coeffRef(snode,mnode)+= penaltyStifness *(Dd2Dxs*Dd2Dym+d2*DDd2DxsDym);
       
         Hiyjx.coeffRef(snode,mnode)+= penaltyStifness *(Dd2Dys*Dd2Dxm+d2*DDd2DysDxm);
         Hiyjy.coeffRef(snode,mnode)+= penaltyStifness *(Dd2Dys*Dd2Dym+d2*DDd2DysDym);
      
         Hixjx.coeffRef(mnode,snode)+= penaltyStifness *(Dd2Dxm*Dd2Dxs+d2*DDd2DxsDxm);
         Hixjy.coeffRef(mnode,snode)+= penaltyStifness *(Dd2Dxm*Dd2Dys+d2*DDd2DxsDym);

         Hiyjx.coeffRef(mnode,snode)+= penaltyStifness *(Dd2Dym*Dd2Dxs+d2*DDd2DxsDym);
         Hiyjy.coeffRef(mnode,snode)+= penaltyStifness *(Dd2Dym*Dd2Dys+d2*DDd2DysDym);
     
     }else if (snode < baseData.numOriginalNodes || mnode < baseData.numOriginalNodes){
         
         
         //replace the ghost node with the original node number
           if (mnode >= baseData.numOriginalNodes) {
               mnode = mnode % baseData.numOriginalNodes;
           }else{
               snode = snode % baseData.numOriginalNodes;
           }
           
           // And then add 0.5 of the Hessian value assuming that this triplet of nodes will be visited again from the other side. If there are cases where this assuption does not hold then the Hessian will be incorrect.
     
         Hixjx.coeffRef(snode,snode)+= penaltyStifness *(Dd2Dxs*Dd2Dxs+d2*DDd2DxsDxs) * 0.5;
         Hixjy.coeffRef(snode,snode)+= penaltyStifness *(Dd2Dxs*Dd2Dys+d2*DDd2DxsDys);
         Hiyjx.coeffRef(snode,snode)+= penaltyStifness *(Dd2Dys*Dd2Dxs+d2*DDd2DxsDys) * 0.5;
         Hiyjy.coeffRef(snode,snode)+= penaltyStifness *(Dd2Dys*Dd2Dys+d2*DDd2DysDys) * 0.5;
         
         Hixjx.coeffRef(mnode,mnode)+= penaltyStifness *(Dd2Dxm*Dd2Dxm+d2*DDd2DxmDxm) * 0.5;
         Hixjy.coeffRef(mnode,mnode)+= penaltyStifness *(Dd2Dxm*Dd2Dym+d2*DDd2DxmDym) * 0.5;
         Hiyjx.coeffRef(mnode,mnode)+= penaltyStifness *(Dd2Dym*Dd2Dxm+d2*DDd2DxmDym) * 0.5;
         Hiyjy.coeffRef(mnode,mnode)+= penaltyStifness *(Dd2Dym*Dd2Dym+d2*DDd2DymDym) * 0.5;
         
     }
    
    
    
    
    
    //____________________________________________________this was the original version of the code where I ignored ghost nodes.
    
//    if ( snode < baseData.numOriginalNodes && mnode < baseData.numOriginalNodes){
//        Hixjx.coeffRef(snode,snode)+= -penaltyStifness *(Dd2Dxs*Dd2Dxs+d2*DDd2DxsDxs);
//        Hixjy.coeffRef(snode,snode)+= -penaltyStifness *(Dd2Dxs*Dd2Dys+d2*DDd2DxsDys);
//        Hiyjx.coeffRef(snode,snode)+= -penaltyStifness *(Dd2Dys*Dd2Dxs+d2*DDd2DxsDys);
//        Hiyjy.coeffRef(snode,snode)+= -penaltyStifness *(Dd2Dys*Dd2Dys+d2*DDd2DysDys);
//    }
//    
//    if ( mnode < baseData.numOriginalNodes){
//        Hixjx.coeffRef(mnode,mnode)+= -penaltyStifness *(Dd2Dxm*Dd2Dxm+d2*DDd2DxmDxm);
//        Hixjy.coeffRef(mnode,mnode)+= -penaltyStifness *(Dd2Dxm*Dd2Dym+d2*DDd2DxmDym);
//        Hiyjx.coeffRef(mnode,mnode)+= -penaltyStifness *(Dd2Dym*Dd2Dxm+d2*DDd2DxmDym);
//        Hiyjy.coeffRef(mnode,mnode)+= -penaltyStifness *(Dd2Dym*Dd2Dym+d2*DDd2DymDym);
//    }
       
}
