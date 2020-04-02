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
#include <stdio.h>

void Configuration::add_d1_contributions_to_Hessian(double penaltyStifness, double xs,double ys,double x0,double y0,double x1,double y1, int snode, int node0, int node1, const BaseSysData& baseData){

double L = pow((pow((x1 - x0),2) + pow((y1 - y0),2)),0.5);
double d1 = (xs - x0)*(y0 - y1)/L + (ys - y0)*(x1 - x0)/L;  // this is the equation used in te derivations. Its the negative of the signedGap equation in NTS function
 
double Dd1Dxs=(y0 - y1)/pow(pow(-x0 + x1,2) + pow(-y0 + y1,2),0.5);

double Dd1Dys=(-x0 + x1)/pow(pow(-x0 + x1,2) + pow(-y0 + y1,2),0.5);

double Dd1Dx0=(1.*(x0 - 1.*x1)*(x0 - 1.*xs)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (-y0 + y1)/pow(pow(x0 - x1,2) + pow(y0 - y1,2),0.5) + (y0 - ys)/pow(pow(x0 - x1,2) + pow(y0 - y1,2),0.5) - (1.*pow(x0 - 1.*x1,2)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double Dd1Dy0=(x0 - x1)/pow(pow(x0 - x1,2) + pow(y0 - y1,2),0.5) + (-x0 + xs)/pow(pow(x0 - x1,2) + pow(y0 - y1,2),0.5) + (1.*(x0 - 1.*xs)*pow(y0 - 1.*y1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) - (1.*(x0 - 1.*x1)*(y0 - 1.*y1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double Dd1Dx1=(-1.*(x0 - 1.*x1)*(x0 - 1.*xs)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (1.*pow(x0 - 1.*x1,2)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (-y0 + ys)/pow(pow(x0 - x1,2) + pow(y0 - y1,2),0.5);

double Dd1Dy1=(x0 - xs)/pow(pow(x0 - x1,2) + pow(y0 - y1,2),0.5) - (1.*(x0 - 1.*xs)*pow(y0 - 1.*y1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (1.*(x0 - 1.*x1)*(y0 - 1.*y1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double DDd1DxsDxs=0;

double DDd1DxsDys=0;

double DDd1DxsDx0=(-1.*(x0 - 1.*x1)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double DDd1DxsDy0=pow(pow(x0 - x1,2) + pow(y0 - y1,2),-0.5) - (1.*pow(y0 - 1.*y1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double DDd1DxsDx1=(1.*(x0 - 1.*x1)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double DDd1DxsDy1=-pow(pow(x0 - x1,2) + pow(y0 - y1,2),-0.5) + (1.*pow(y0 - 1.*y1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double DDd1DysDys=0;

double DDd1DysDx0=(1.*pow(x0 - 1.*x1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) - pow(pow(x0 - x1,2) + pow(y0 - y1,2),-0.5);

double DDd1DysDy0=(1.*(x0 - 1.*x1)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double DDd1DysDx1=(-1.*pow(x0 - 1.*x1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + pow(pow(x0 - x1,2) + pow(y0 - y1,2),-0.5);

double DDd1DysDy1=(-1.*(x0 - 1.*x1)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double DDd1Dx0Dx0=(-3.*pow(x0 - 1.*x1,2)*(x0 - 1.*xs)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) + (2.*(x0 - 1.*x1)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (1.*(x0 - 1.*xs)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (3.*pow(x0 - 1.*x1,3)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) - (3.*(x0 - 1.*x1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double DDd1Dx0Dy0=(-1.*pow(x0 - 1.*x1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (1.*(x0 - 1.*x1)*(x0 - 1.*xs))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) - (3.*(x0 - 1.*x1)*(x0 - 1.*xs)*pow(y0 - 1.*y1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) + (1.*pow(y0 - 1.*y1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (3.*pow(x0 - 1.*x1,2)*(y0 - 1.*y1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) - (1.*(y0 - 1.*y1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double DDd1Dx0Dx1=(3.*pow(x0 - 1.*x1,2)*(x0 - 1.*xs)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) - (1.*(x0 - 1.*x1)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) - (1.*(x0 - 1.*xs)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) - (3.*pow(x0 - 1.*x1,3)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) + (3.*(x0 - 1.*x1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double DDd1Dx0Dy1=(-1.*(x0 - 1.*x1)*(x0 - 1.*xs))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + pow(pow(x0 - x1,2) + pow(y0 - y1,2),-0.5) + (3.*(x0 - 1.*x1)*(x0 - 1.*xs)*pow(y0 - 1.*y1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) - (1.*pow(y0 - 1.*y1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) - (3.*pow(x0 - 1.*x1,2)*(y0 - 1.*y1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) + (1.*(y0 - 1.*y1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double DDd1Dy0Dy0=(-2.*(x0 - 1.*x1)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (3.*(x0 - 1.*xs)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) - (3.*(x0 - 1.*xs)*pow(y0 - 1.*y1,3))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) - (1.*(x0 - 1.*x1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (3.*(x0 - 1.*x1)*pow(y0 - 1.*y1,2)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5);

double DDd1Dy0Dx1=(1.*pow(x0 - 1.*x1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) - (1.*(x0 - 1.*x1)*(x0 - 1.*xs))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) - pow(pow(x0 - x1,2) + pow(y0 - y1,2),-0.5) + (3.*(x0 - 1.*x1)*(x0 - 1.*xs)*pow(y0 - 1.*y1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) - (3.*pow(x0 - 1.*x1,2)*(y0 - 1.*y1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) + (1.*(y0 - 1.*y1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double DDd1Dy0Dy1=(1.*(x0 - 1.*x1)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) - (3.*(x0 - 1.*xs)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (3.*(x0 - 1.*xs)*pow(y0 - 1.*y1,3))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) + (1.*(x0 - 1.*x1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) - (3.*(x0 - 1.*x1)*pow(y0 - 1.*y1,2)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5);

double DDd1Dx1Dx1=(-3.*pow(x0 - 1.*x1,2)*(x0 - 1.*xs)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) + (1.*(x0 - 1.*xs)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (3.*pow(x0 - 1.*x1,3)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) - (3.*(x0 - 1.*x1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double DDd1Dx1Dy1=(1.*(x0 - 1.*x1)*(x0 - 1.*xs))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) - (3.*(x0 - 1.*x1)*(x0 - 1.*xs)*pow(y0 - 1.*y1,2))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) + (3.*pow(x0 - 1.*x1,2)*(y0 - 1.*y1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) - (1.*(y0 - 1.*y1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5);

double DDd1Dy1Dy1=(3.*(x0 - 1.*xs)*(y0 - 1.*y1))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) - (3.*(x0 - 1.*xs)*pow(y0 - 1.*y1,3))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5) - (1.*(x0 - 1.*x1)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),1.5) + (3.*(x0 - 1.*x1)*pow(y0 - 1.*y1,2)*(y0 - 1.*ys))/pow(pow(x0 - x1,2) + pow(y0 - y1,2),2.5);
//    penaltyStifness *= 100;
    
   
    
    // if all nodes are original
   if ( (snode < baseData.numOriginalNodes) &&  (node0 < baseData.numOriginalNodes) ){
       
       Hixjx.coeffRef(snode,snode)+= penaltyStifness *(Dd1Dxs*Dd1Dxs+d1*DDd1DxsDxs);
       Hixjy.coeffRef(snode,snode)+= penaltyStifness *(Dd1Dxs*Dd1Dys+d1*DDd1DxsDys);
       
       Hiyjx.coeffRef(snode,snode)+= penaltyStifness *(Dd1Dys*Dd1Dxs+d1*DDd1DxsDys);
       Hiyjy.coeffRef(snode,snode)+= penaltyStifness *(Dd1Dys*Dd1Dys+d1*DDd1DysDys);
       
       Hixjx.coeffRef(node0,node0)+= penaltyStifness *(Dd1Dx0*Dd1Dx0+d1*DDd1Dx0Dx0);
       Hixjy.coeffRef(node0,node0)+= penaltyStifness *(Dd1Dx0*Dd1Dy0+d1*DDd1Dx0Dy0);
       
       Hiyjx.coeffRef(node0,node0)+= penaltyStifness *(Dd1Dy0*Dd1Dx0+d1*DDd1Dx0Dy0);
       Hiyjy.coeffRef(node0,node0)+= penaltyStifness *(Dd1Dy0*Dd1Dy0+d1*DDd1Dy0Dy0);
       
       Hixjx.coeffRef(node1,node1)+= penaltyStifness *(Dd1Dx1*Dd1Dx1+d1*DDd1Dx1Dx1);
       Hixjy.coeffRef(node1,node1)+= penaltyStifness *(Dd1Dx1*Dd1Dy1+d1*DDd1Dx1Dy1);
       
       Hiyjx.coeffRef(node1,node1)+= penaltyStifness *(Dd1Dy1*Dd1Dx1+d1*DDd1Dx1Dy1);
       Hiyjy.coeffRef(node1,node1)+= penaltyStifness *(Dd1Dy1*Dd1Dy1+d1*DDd1Dy1Dy1);
       
       Hixjx.coeffRef(snode,node0)+= penaltyStifness *(Dd1Dxs*Dd1Dx0+d1*DDd1DxsDx0);
       Hixjy.coeffRef(snode,node0)+= penaltyStifness *(Dd1Dxs*Dd1Dy0+d1*DDd1DxsDy0);
       
       Hiyjx.coeffRef(snode,node0)+= penaltyStifness *(Dd1Dys*Dd1Dx0+d1*DDd1DysDx0);
       Hiyjy.coeffRef(snode,node0)+= penaltyStifness *(Dd1Dys*Dd1Dy0+d1*DDd1DysDy0);
       
       Hixjx.coeffRef(node0,snode)+= penaltyStifness *(Dd1Dx0*Dd1Dxs+d1*DDd1DxsDx0);
       Hixjy.coeffRef(node0,snode)+= penaltyStifness *(Dd1Dx0*Dd1Dys+d1*DDd1DysDx0);
       
       Hiyjx.coeffRef(node0,snode)+= penaltyStifness *(Dd1Dy0*Dd1Dxs+d1*DDd1DxsDy0);
       Hiyjy.coeffRef(node0,snode)+= penaltyStifness *(Dd1Dy0*Dd1Dys+d1*DDd1DysDy0);
       
       Hixjx.coeffRef(snode,node1)+= penaltyStifness *(Dd1Dxs*Dd1Dx1+d1*DDd1DxsDx1);
       Hixjy.coeffRef(snode,node1)+= penaltyStifness *(Dd1Dxs*Dd1Dy1+d1*DDd1DxsDy1);

       Hiyjx.coeffRef(snode,node1)+= penaltyStifness *(Dd1Dys*Dd1Dx1+d1*DDd1DysDx1);
       Hiyjy.coeffRef(snode,node1)+= penaltyStifness *(Dd1Dys*Dd1Dy1+d1*DDd1DysDy1);
       
       Hixjx.coeffRef(node1,snode)+= penaltyStifness *(Dd1Dx1*Dd1Dxs+d1*DDd1DxsDx1);
       Hixjy.coeffRef(node1,snode)+= penaltyStifness *(Dd1Dx1*Dd1Dys+d1*DDd1DysDx1);
       
       Hiyjx.coeffRef(node1,snode)+= penaltyStifness *(Dd1Dy1*Dd1Dxs+d1*DDd1DxsDy1);
       Hiyjy.coeffRef(node1,snode)+= penaltyStifness *(Dd1Dy1*Dd1Dys+d1*DDd1DysDy1);
       
       
       Hixjx.coeffRef(node0,node1)+= penaltyStifness *(Dd1Dx0*Dd1Dx1+d1*DDd1Dx0Dx1);
       Hixjy.coeffRef(node0,node1)+= penaltyStifness *(Dd1Dx0*Dd1Dy1+d1*DDd1Dx0Dy1);
       
       Hiyjx.coeffRef(node0,node1)+= penaltyStifness *(Dd1Dy0*Dd1Dx1+d1*DDd1Dy0Dx1);
       Hiyjy.coeffRef(node0,node1)+= penaltyStifness *(Dd1Dy0*Dd1Dy1+d1*DDd1Dy0Dy1);
       
       Hixjx.coeffRef(node1,node0)+= penaltyStifness *(Dd1Dx1*Dd1Dx0+d1*DDd1Dx0Dx1);
       Hixjy.coeffRef(node1,node0)+= penaltyStifness *(Dd1Dx1*Dd1Dy0+d1*DDd1Dy0Dx1);
       
       Hiyjx.coeffRef(node1,node0)+= penaltyStifness *(Dd1Dy1*Dd1Dx0+d1*DDd1Dx0Dy1);
       Hiyjy.coeffRef(node1,node0)+= penaltyStifness *(Dd1Dy1*Dd1Dy0+d1*DDd1Dy0Dy1);
       
       
   }else if (  snode < baseData.numOriginalNodes ||  node0 < baseData.numOriginalNodes  ){
       
       
       //replace the ghost node with the original node number
       if (node0 >= baseData.numOriginalNodes) {
           node0 = node0 % baseData.numOriginalNodes;
           node1 = node1 % baseData.numOriginalNodes;
       }else{
           snode = snode % baseData.numOriginalNodes;
       }
       
       // And then add 0.5 of the Hessian value assuming that this triplet of nodes will be visited again from the other side. If there are cases where this assuption does not hold then the Hessian will be incorrect.
       
       Hixjx.coeffRef(snode,snode)+= penaltyStifness *(Dd1Dxs*Dd1Dxs+d1*DDd1DxsDxs) * 0.5;
       Hixjy.coeffRef(snode,snode)+= penaltyStifness *(Dd1Dxs*Dd1Dys+d1*DDd1DxsDys) * 0.5;
       
       Hiyjx.coeffRef(snode,snode)+= penaltyStifness *(Dd1Dys*Dd1Dxs+d1*DDd1DxsDys) * 0.5;
       Hiyjy.coeffRef(snode,snode)+= penaltyStifness *(Dd1Dys*Dd1Dys+d1*DDd1DysDys) * 0.5;
       
       Hixjx.coeffRef(node0,node0)+= penaltyStifness *(Dd1Dx0*Dd1Dx0+d1*DDd1Dx0Dx0) * 0.5;
       Hixjy.coeffRef(node0,node0)+= penaltyStifness *(Dd1Dx0*Dd1Dy0+d1*DDd1Dx0Dy0) * 0.5;
       
       Hiyjx.coeffRef(node0,node0)+= penaltyStifness *(Dd1Dy0*Dd1Dx0+d1*DDd1Dx0Dy0) * 0.5;
       Hiyjy.coeffRef(node0,node0)+= penaltyStifness *(Dd1Dy0*Dd1Dy0+d1*DDd1Dy0Dy0) * 0.5;
       
       Hixjx.coeffRef(node1,node1)+= penaltyStifness *(Dd1Dx1*Dd1Dx1+d1*DDd1Dx1Dx1) * 0.5;
       Hixjy.coeffRef(node1,node1)+= penaltyStifness *(Dd1Dx1*Dd1Dy1+d1*DDd1Dx1Dy1) * 0.5;
       
       Hiyjx.coeffRef(node1,node1)+= penaltyStifness *(Dd1Dy1*Dd1Dx1+d1*DDd1Dx1Dy1) * 0.5;
       Hiyjy.coeffRef(node1,node1)+= penaltyStifness *(Dd1Dy1*Dd1Dy1+d1*DDd1Dy1Dy1) * 0.5;
       
       Hixjx.coeffRef(snode,node0)+= penaltyStifness *(Dd1Dxs*Dd1Dx0+d1*DDd1DxsDx0) * 0.5;
       Hixjy.coeffRef(snode,node0)+= penaltyStifness *(Dd1Dxs*Dd1Dy0+d1*DDd1DxsDy0) * 0.5;
       
       Hiyjx.coeffRef(snode,node0)+= penaltyStifness *(Dd1Dys*Dd1Dx0+d1*DDd1DysDx0) * 0.5;
       Hiyjy.coeffRef(snode,node0)+= penaltyStifness *(Dd1Dys*Dd1Dy0+d1*DDd1DysDy0) * 0.5;
       
       Hixjx.coeffRef(node0,snode)+= penaltyStifness *(Dd1Dx0*Dd1Dxs+d1*DDd1DxsDx0) * 0.5;
       Hixjy.coeffRef(node0,snode)+= penaltyStifness *(Dd1Dx0*Dd1Dys+d1*DDd1DysDx0) * 0.5;
       
       Hiyjx.coeffRef(node0,snode)+= penaltyStifness *(Dd1Dy0*Dd1Dxs+d1*DDd1DxsDy0) * 0.5;
       Hiyjy.coeffRef(node0,snode)+= penaltyStifness *(Dd1Dy0*Dd1Dys+d1*DDd1DysDy0) * 0.5;
       
       Hixjx.coeffRef(snode,node1)+= penaltyStifness *(Dd1Dxs*Dd1Dx1+d1*DDd1DxsDx1) * 0.5;
       Hixjy.coeffRef(snode,node1)+= penaltyStifness *(Dd1Dxs*Dd1Dy1+d1*DDd1DxsDy1) * 0.5;

       Hiyjx.coeffRef(snode,node1)+= penaltyStifness *(Dd1Dys*Dd1Dx1+d1*DDd1DysDx1) * 0.5;
       Hiyjy.coeffRef(snode,node1)+= penaltyStifness *(Dd1Dys*Dd1Dy1+d1*DDd1DysDy1) * 0.5;
       
       Hixjx.coeffRef(node1,snode)+= penaltyStifness *(Dd1Dx1*Dd1Dxs+d1*DDd1DxsDx1) * 0.5;
       Hixjy.coeffRef(node1,snode)+= penaltyStifness *(Dd1Dx1*Dd1Dys+d1*DDd1DysDx1) * 0.5;
       
       Hiyjx.coeffRef(node1,snode)+= penaltyStifness *(Dd1Dy1*Dd1Dxs+d1*DDd1DxsDy1) * 0.5;
       Hiyjy.coeffRef(node1,snode)+= penaltyStifness *(Dd1Dy1*Dd1Dys+d1*DDd1DysDy1) * 0.5;
       
       
       Hixjx.coeffRef(node0,node1)+= penaltyStifness *(Dd1Dx0*Dd1Dx1+d1*DDd1Dx0Dx1) * 0.5;
       Hixjy.coeffRef(node0,node1)+= penaltyStifness *(Dd1Dx0*Dd1Dy1+d1*DDd1Dx0Dy1) * 0.5;
       
       Hiyjx.coeffRef(node0,node1)+= penaltyStifness *(Dd1Dy0*Dd1Dx1+d1*DDd1Dy0Dx1) * 0.5;
       Hiyjy.coeffRef(node0,node1)+= penaltyStifness *(Dd1Dy0*Dd1Dy1+d1*DDd1Dy0Dy1) * 0.5;
       
       Hixjx.coeffRef(node1,node0)+= penaltyStifness *(Dd1Dx1*Dd1Dx0+d1*DDd1Dx0Dx1) * 0.5;
       Hixjy.coeffRef(node1,node0)+= penaltyStifness *(Dd1Dx1*Dd1Dy0+d1*DDd1Dy0Dx1) * 0.5;
       
       Hiyjx.coeffRef(node1,node0)+= penaltyStifness *(Dd1Dy1*Dd1Dx0+d1*DDd1Dx0Dy1) * 0.5;
       Hiyjy.coeffRef(node1,node0)+= penaltyStifness *(Dd1Dy1*Dd1Dy0+d1*DDd1Dy0Dy1) * 0.5;
   }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   //____________________________________________________this was the original version of the code where I ignored ghost nodes.
//    if ( snode < baseData.numOriginalNodes){
//
//
//
//        Hixjx.coeffRef(snode,snode)+= -penaltyStifness *(Dd1Dxs*Dd1Dxs+d1*DDd1DxsDxs);
//        Hixjy.coeffRef(snode,snode)+= -penaltyStifness *(Dd1Dxs*Dd1Dys+d1*DDd1DxsDys);
//
//        Hiyjx.coeffRef(snode,snode)+= -penaltyStifness *(Dd1Dys*Dd1Dxs+d1*DDd1DxsDys);
//        Hiyjy.coeffRef(snode,snode)+= -penaltyStifness *(Dd1Dys*Dd1Dys+d1*DDd1DysDys);
//    }
//
//    if ( node0 < baseData.numOriginalNodes){
//
//        Hixjx.coeffRef(node0,node0)+= -penaltyStifness *(Dd1Dx0*Dd1Dx0+d1*DDd1Dx0Dx0);
//        Hixjy.coeffRef(node0,node0)+= -penaltyStifness *(Dd1Dx0*Dd1Dy0+d1*DDd1Dx0Dy0);
//
//        Hiyjx.coeffRef(node0,node0)+= -penaltyStifness *(Dd1Dy0*Dd1Dx0+d1*DDd1Dx0Dy0);
//        Hiyjy.coeffRef(node0,node0)+= -penaltyStifness *(Dd1Dy0*Dd1Dy0+d1*DDd1Dy0Dy0);
//    }
//
//    if ( node1 < baseData.numOriginalNodes){
//
//        Hixjx.coeffRef(node1,node1)+= -penaltyStifness *(Dd1Dx1*Dd1Dx1+d1*DDd1Dx1Dx1);
//        Hixjy.coeffRef(node1,node1)+= -penaltyStifness *(Dd1Dx1*Dd1Dy1+d1*DDd1Dx1Dy1);
//
//        Hiyjx.coeffRef(node1,node1)+= -penaltyStifness *(Dd1Dy1*Dd1Dx1+d1*DDd1Dx1Dy1);
//        Hiyjy.coeffRef(node1,node1)+= -penaltyStifness *(Dd1Dy1*Dd1Dy1+d1*DDd1Dy1Dy1);
//    }
    
//     if ( snode < baseData.numOriginalNodes && node0 < baseData.numOriginalNodes){
//
//         Hixjx.coeffRef(snode,node0)+= -penaltyStifness *(Dd1Dxs*Dd1Dx0+d1*DDd1DxsDx0);
//         Hixjy.coeffRef(snode,node0)+= -penaltyStifness *(Dd1Dxs*Dd1Dy0+d1*DDd1DxsDy0);
//
//         Hiyjx.coeffRef(snode,node0)+= -penaltyStifness *(Dd1Dys*Dd1Dx0+d1*DDd1DysDx0);
//         Hiyjy.coeffRef(snode,node0)+= -penaltyStifness *(Dd1Dys*Dd1Dy0+d1*DDd1DysDy0);
//
//         Hixjx.coeffRef(node0,snode)+= -penaltyStifness *(Dd1Dx0*Dd1Dxs+d1*DDd1DxsDx0);
//         Hixjy.coeffRef(node0,snode)+= -penaltyStifness *(Dd1Dx0*Dd1Dys+d1*DDd1DysDx0);
//
//         Hiyjx.coeffRef(node0,snode)+= -penaltyStifness *(Dd1Dy0*Dd1Dxs+d1*DDd1DxsDy0);
//         Hiyjy.coeffRef(node0,snode)+= -penaltyStifness *(Dd1Dy0*Dd1Dys+d1*DDd1DysDy0);
//     }
//
//    if ( snode < baseData.numOriginalNodes && node1 < baseData.numOriginalNodes){
//
//         Hixjx.coeffRef(snode,node1)+= -penaltyStifness *(Dd1Dxs*Dd1Dx1+d1*DDd1DxsDx1);
//         Hixjy.coeffRef(snode,node1)+= -penaltyStifness *(Dd1Dxs*Dd1Dy1+d1*DDd1DxsDy1);
//
//         Hiyjx.coeffRef(snode,node1)+= -penaltyStifness *(Dd1Dys*Dd1Dx1+d1*DDd1DysDx1);
//         Hiyjy.coeffRef(snode,node1)+= -penaltyStifness *(Dd1Dys*Dd1Dy1+d1*DDd1DysDy1);
//
//         Hixjx.coeffRef(node1,snode)+= -penaltyStifness *(Dd1Dx1*Dd1Dxs+d1*DDd1DxsDx1);
//         Hixjy.coeffRef(node1,snode)+= -penaltyStifness *(Dd1Dx1*Dd1Dys+d1*DDd1DysDx1);
//
//         Hiyjx.coeffRef(node1,snode)+= -penaltyStifness *(Dd1Dy1*Dd1Dxs+d1*DDd1DxsDy1);
//         Hiyjy.coeffRef(node1,snode)+= -penaltyStifness *(Dd1Dy1*Dd1Dys+d1*DDd1DysDy1);
//    }
    
    
    
//     if ( node0 < baseData.numOriginalNodes && node1 < baseData.numOriginalNodes){
//
//          Hixjx.coeffRef(node0,node1)+= -penaltyStifness *(Dd1Dx0*Dd1Dx1+d1*DDd1Dx0Dx1);
//          Hixjy.coeffRef(node0,node1)+= -penaltyStifness *(Dd1Dx0*Dd1Dy1+d1*DDd1Dx0Dy1);
//
//          Hiyjx.coeffRef(node0,node1)+= -penaltyStifness *(Dd1Dy0*Dd1Dx1+d1*DDd1Dy0Dx1);
//          Hiyjy.coeffRef(node0,node1)+= -penaltyStifness *(Dd1Dy0*Dd1Dy1+d1*DDd1Dy0Dy1);
//
//          Hixjx.coeffRef(node1,node0)+= -penaltyStifness *(Dd1Dx1*Dd1Dx0+d1*DDd1Dx0Dx1);
//          Hixjy.coeffRef(node1,node0)+= -penaltyStifness *(Dd1Dx1*Dd1Dy0+d1*DDd1Dy0Dx1);
//
//          Hiyjx.coeffRef(node1,node0)+= -penaltyStifness *(Dd1Dy1*Dd1Dx0+d1*DDd1Dx0Dy1);
//          Hiyjy.coeffRef(node1,node0)+= -penaltyStifness *(Dd1Dy1*Dd1Dy0+d1*DDd1Dy0Dy1);
//
//     }
    
      
}

