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


void Configuration::calculate_the_Hessian_H(const Parameters &pars){
    
    // First construct the stiffness matrix
    calculate_monolithic_stiffness_tensor_K(pars);
    
    // Now so the first tensor multiplication
    for (int colID=0; colID<gradX.cols(); colID++) {
       
        KMxxjx.col(colID) = gradX.col(colID).cwiseProduct(Kxxxx) + gradY.col(colID).cwiseProduct(Kxxxy);
        KMxxjy.col(colID) = gradX.col(colID).cwiseProduct(Kxxyx) + gradY.col(colID).cwiseProduct(Kxxyy);
        KMxyjx.col(colID) = gradX.col(colID).cwiseProduct(Kxxxy) + gradY.col(colID).cwiseProduct(Kxyxy);
        KMxyjy.col(colID) = gradX.col(colID).cwiseProduct(Kxyyx) + gradY.col(colID).cwiseProduct(Kxyyy);
        KMyxjx.col(colID) = gradX.col(colID).cwiseProduct(Kxxyx) + gradY.col(colID).cwiseProduct(Kxyyx);
        KMyxjy.col(colID) = gradX.col(colID).cwiseProduct(Kyxyx) + gradY.col(colID).cwiseProduct(Kyxyy);
        KMyyjx.col(colID) = gradX.col(colID).cwiseProduct(Kxxyy) + gradY.col(colID).cwiseProduct(Kxyyy);
        KMyyjy.col(colID) = gradX.col(colID).cwiseProduct(Kyxyy) + gradY.col(colID).cwiseProduct(Kyyyy);
        
        // multiply by the refArea of each element
        KMxxjx.col(colID) = KMxxjx.col(colID).cwiseProduct(refArea);
        KMxxjy.col(colID) = KMxxjy.col(colID).cwiseProduct(refArea);
        KMxyjx.col(colID) = KMxyjx.col(colID).cwiseProduct(refArea);
        KMxyjy.col(colID) = KMxyjy.col(colID).cwiseProduct(refArea);
        KMyxjx.col(colID) = KMyxjx.col(colID).cwiseProduct(refArea);
        KMyxjy.col(colID) = KMyxjy.col(colID).cwiseProduct(refArea);
        KMyyjx.col(colID) = KMyyjx.col(colID).cwiseProduct(refArea);
        KMyyjy.col(colID) = KMyyjy.col(colID).cwiseProduct(refArea);
    }
    
    
    
    //DCalculate initial H (that is before multipying it by refareas)
    Hixjx = gradX.transpose()*KMxxjx + gradY.transpose()*KMxyjx;
    Hixjy = gradX.transpose()*KMxxjy + gradY.transpose()*KMxyjy;
    Hiyjx = gradX.transpose()*KMyxjx + gradY.transpose()*KMyyjx;
    Hiyjy = gradX.transpose()*KMyxjy + gradY.transpose()*KMyyjy;
    
   
}
