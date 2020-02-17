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


void Configuration::calculate_monolithic_stiffness_tensor_K(const Parameters& pars){
    
    WmPrimePrimeJSquared = pars.kTOverOmega * (  ( 2*pars.chi*(1-areaRatio.array()) + areaRatio.array() ) / areaRatio.array() / (areaRatio.array() - 1) );
    
    prefactorA = pars.NkT + WmPrimePrimeJSquared.array();
    prefactorB = pars.NkT - swellingPressurePerEle.array() * areaRatio.array() ;
    prefactorC = WmPrimePrimeJSquared.array() + swellingPressurePerEle.array() * areaRatio.array() ;
    
    Kxxxx = pars.NkT + prefactorA.array() * invDefGradTransXX.array()* invDefGradTransXX.array();
    Kxxxy = prefactorA.array() * invDefGradTransYX.array()* invDefGradTransXX.array();
    Kxxyx = prefactorA.array() * invDefGradTransXX.array()*invDefGradTransXY.array();
    Kxxyy = prefactorB.array() * invDefGradTransYX.array()*invDefGradTransXY.array() + prefactorC.array() * invDefGradTransYY.array() * invDefGradTransXX.array();
    Kxyxy  = pars.NkT + prefactorA.array() * invDefGradTransYX.array() * invDefGradTransYX .array();
    Kxyyx = prefactorB.array() * invDefGradTransXX.array()*invDefGradTransYY.array() + prefactorC.array()*invDefGradTransXY.array()*invDefGradTransYX.array();
    Kxyyy = prefactorA.array()*invDefGradTransYX.array()*invDefGradTransYY.array();
    Kyxyx = pars.NkT + prefactorA.array()*invDefGradTransXY.array()*invDefGradTransXY.array();
    Kyxyy = prefactorA.array()*invDefGradTransYY.array()*invDefGradTransXY.array();
    Kyyyy = pars.NkT + prefactorA.array()*invDefGradTransYY.array()*invDefGradTransYY.array();
}

