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


struct closestPointOnHermitianCurve find_closest_point_on_Hermitian_interpolation(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double alpha){
    
    struct closestPointOnHermitianCurve soluionPoint;
    double g=-1;
    double xg,yg, DxDg, DyDg, DDxDDg, DDyDDg, F, F_prime;
    int i=0;
    soluionPoint.gap = 999;
   do{
       i++;

        xg = (alpha*(1 - 2*g + pow(g,2))*x2)/4. + ((4 - 2*alpha*(1 + pow(g,2)))*x3)/4. + (alpha*(1 + 2*g + pow(g,2))*x4)/4.;
        yg = (alpha*(1 - 2*g + pow(g,2))*y2)/4. + ((4 - 2*alpha*(1 + pow(g,2)))*y3)/4. + (alpha*(1 + 2*g + pow(g,2))*y4)/4.;
        DxDg = (alpha*(-2 + 2*g)*x2)/4. - alpha*g*x3 + (alpha*(2 + 2*g)*x4)/4.;
        DyDg = (alpha*(-2 + 2*g)*y2)/4. - alpha*g*y3 + (alpha*(2 + 2*g)*y4)/4.;
        DDxDDg = (alpha*x2)/2. - alpha*x3 + (alpha*x4)/2.;
        DDyDDg = (alpha*y2)/2. - alpha*y3 + (alpha*y4)/2.;
        F = (xg-x1)*DxDg+(yg-y1)*DyDg;
        F_prime = (DxDg*DxDg+DyDg*DyDg) + ((xg-x1)*DDxDDg+(yg-y1)*DDyDDg);
        g = g - F/F_prime;


   } while(abs(F) > 1E-12);
    
//    if (i<10000){
    soluionPoint.g= g + F/F_prime;
    soluionPoint.x = xg;
    soluionPoint.y = yg;
    soluionPoint.nx = DyDg/sqrt(DxDg*DxDg+DyDg*DyDg);
    soluionPoint.ny = - DxDg/sqrt(DxDg*DxDg+DyDg*DyDg);
    if (g<=1 && g>=-1) soluionPoint.gap =(x1-xg)*soluionPoint.nx + (y1-yg)*soluionPoint.ny ;
//    }else{
//        soluionPoint.g= 99.;
//        soluionPoint.x = xg;
//        soluionPoint.y = yg;
//        soluionPoint.nx = DyDg/sqrt(DxDg*DxDg+DyDg*DyDg);
//        soluionPoint.ny = - DxDg/sqrt(DxDg*DxDg+DyDg*DyDg);
//    }
    
    
//    std::cout << "it  " << i << std::endl;
//    std::cout << "abs(F)  " << abs(F) << std::endl;
    return soluionPoint;
}

struct point HermitianInterpolation(double x2, double x3, double x4, double y2, double y3, double y4, double alpha, double g){
    
    struct point soluionPoint;
    soluionPoint.x = (alpha*(1 - 2*g + pow(g,2))*x2)/4. + ((4 - 2*alpha*(1 + pow(g,2)))*x3)/4. + (alpha*(1 + 2*g + pow(g,2))*x4)/4.;
    soluionPoint.y = (alpha*(1 - 2*g + pow(g,2))*y2)/4. + ((4 - 2*alpha*(1 + pow(g,2)))*y3)/4. + (alpha*(1 + 2*g + pow(g,2))*y4)/4.;
  
    return soluionPoint;
}


void testFire2(const Parameters& pars){
    
    
    
    // initialize some varialbles
        double FIRE_alpha = pars.FIRE_alpha_start;
        long FIRE_N_positive=0;
        long FIRE_N_negative=0;
        double FIRE_dt = pars.FIRE_dt_start;
        double FdotF;
        double VdotV;
        double scale1=0;
        double scale2=0;
        double x=10,y=10;
        const double Pi=3.14;
        double power, L2R;
        // initialize the system
        double fx =0;
        double fy = 0;
        // initialize velocities
        double vx =0, vy=0;
        
       
        // Leap Frog integration initialization
        if (pars.integrator == 1){
            vx -= 0.5*FIRE_dt *fx;
            vy -= 0.5*FIRE_dt *fy;
        }


        for (long i=1; i<= pars.FIRE_Nmax; i++) {
            

            fx = x/5. + ((Pi*x)/sqrt(pow(x,2) + pow(y,2)) - y/(2.*pow(x,2)*(1 + pow(y,2)/pow(x,2))))*cos(Pi*sqrt(pow(x,2) + pow(y,2)) + atan(y/x)/2.);
            fy = y/5. + ((Pi*y)/sqrt(pow(x,2) + pow(y,2)) + 1/(2.*x*(1 + pow(y,2)/pow(x,2))))*cos(Pi*sqrt(pow(x,2) + pow(y,2)) + atan(y/x)/2.);
            
            L2R = sqrt(fx*fx + fy*fy);
            
            
            
            if (L2R <= pars.FIRE_RTolerance){
                std::cout << " Done !" << std::endl;
                break;
            }else if ( isnan(L2R) || isnan(fx) ||  isnan(fy)){
                std::cout << "System blew up !" << std::endl;
                exit(1);
            }
            
             power = fx*vx + fy*vy;

                std::cout << "step  " << i << "\n" << std::endl;
                std::cout << "energy  " <<  (pow(x,2) + pow(y,2))/10. + sin(Pi*sqrt(pow(x,2) + pow(y,2)) + atan(y/x)/2.) <<std::endl;
                std::cout << "power   " << power <<std::endl;
                std::cout << "L2R  " << L2R << std::endl;
                std::cout << "FIRE_dt   " << FIRE_dt <<std::endl;
                std::cout << "FIRE_alpha   " << FIRE_alpha  <<std::endl;
                std::cout << "FIRE_Np+   " << FIRE_N_positive  <<std::endl;
                std::cout << "FIRE_Np-   " << FIRE_N_negative  <<std::endl;
                std::cout << "\n" << std::endl;
            
            
            if (power >0) {
                FIRE_N_positive +=1;
                FIRE_N_negative =0;
                if (FIRE_N_positive > pars.FIRE_N_positive_min){
                    FIRE_dt = fmin(FIRE_dt * pars.FIRE_finc, pars.FIRE_dtmax);
                    FIRE_alpha *= pars.FIRE_falpha;
                }
                scale1=(1-FIRE_alpha);
                FdotF = fx*fx+fy*fy;
                VdotV = vx*vx+vy*vy;
                if (FdotF <= 1e-20) scale2 = 0.0;
                else scale2 = FIRE_alpha * sqrt(VdotV/FdotF);
                
            } else {
                
                FIRE_N_positive =0;
                FIRE_N_negative +=1;
                
                if (FIRE_N_negative > pars.FIRE_N_negative_max){
                    std::cout << " Failed to converge !" << std::endl;
                    exit(1);
                }
                
                if (!(pars.FIRE_intialdelay && i < pars.FIRE_N_positive_min)){
                    if ( FIRE_dt*pars.FIRE_fdec >= pars.FIRE_dtmin){
                        FIRE_dt*=pars.FIRE_fdec;
                    }
                    FIRE_alpha = pars.FIRE_alpha_start;
                }
                
                x-= 0.5*vx*FIRE_dt;
                y -= 0.5* vx*FIRE_dt;
                vx=0;
                vy=0;
                
            }
//            vmax= sqrt((mainSys.velocityX.array()*mainSys.velocityX.array()+mainSys.velocityY.array()*mainSys.velocityY.array()).maxCoeff());
//
//            if (vmax*FIRE_dt>pars.verletCellCutoff*0.5){
//                FIRE_dt = pars.verletCellCutoff*0.5/vmax;
//            }
            //MD integration
            if (pars.integrator==0) {
                vx+= fx* FIRE_dt;
                vy+= fy* FIRE_dt;
                
            
                if (power>0.0){
                    vx = scale1 *vx+ scale2 * fx;
                    vy = scale1 * vy+ scale2 *fy;
                }
                x += vx* FIRE_dt;
                y += vy* FIRE_dt;
            }
          
            

            
        }
            
    }






////////////////////

void testFire(const Parameters& pars){
    
    
    
    // initialize some varialbles
        double FIRE_alpha = pars.FIRE_alpha_start;
        long FIRE_N_positive=0;
        double FIRE_dt = pars.FIRE_dt_start;
        double FdotF;
        double VdotV;
        double scale1=0;
        double scale2=0;
        double x=10,y=0;
        const double Pi=3.14;
        double power, L2R;
        // initialize the system
        double fx =0;
        double fy = 0;
        // initialize velocities
        double vx =0, vy=0;
        
       
        // Leap Frog integration initialization
        if (pars.integrator == 1){
            vx -= 0.5*FIRE_dt *fx;
            vy -= 0.5*FIRE_dt *fy;
        }


        for (long i=1; i<= pars.FIRE_Nmax; i++) {
            

//            fx = x/5. + ((Pi*x)/sqrt(pow(x,2) + pow(y,2)) - y/(2.*pow(x,2)*(1 + pow(y,2)/pow(x,2))))*cos(Pi*sqrt(pow(x,2) + pow(y,2)) + atan(y/x)/2.);
//            fy = y/5. + ((Pi*y)/sqrt(pow(x,2) + pow(y,2)) + 1/(2.*x*(1 + pow(y,2)/pow(x,2))))*cos(Pi*sqrt(pow(x,2) + pow(y,2)) + atan(y/x)/2.);
            
            fx = x/5. + Pi*cos(Pi*x + y/2.);
            fy = cos(Pi*y + y/2.)/2.;

            vx+= fx* FIRE_dt;
            vy+= fy* FIRE_dt;
            
            x += vx* FIRE_dt;
            y += vy* FIRE_dt;
            
            L2R = sqrt(fx*fx + fy*fy);
            
            
            
            if (L2R <= pars.FIRE_RTolerance){
                std::cout << " Done !" << std::endl;
                break;
            }else if ( isnan(L2R) || isnan(fx) ||  isnan(fy)){
                std::cout << "System blew up !" << std::endl;
                exit(1);
            }
            
            power = fx*vx + fy*vy;
            
            FdotF = fx*fx+fy*fy;
            VdotV = vx*vx+vy*vy;
            scale1=(1-FIRE_alpha);
            scale2 = FIRE_alpha * sqrt(VdotV/FdotF);
            vx = scale1 *vx+ scale2 * fx;
            vy = scale1 * vy+ scale2 *fy;

                std::cout << "step  " << i << "\n" << std::endl;
                std::cout << "energy  " << pow(x,2)/10. + sin(Pi*x + y/2.) <<std::endl;
                std::cout << "power   " << power <<std::endl;
                std::cout << "L2R  " << L2R << std::endl;
//                std::cout << "FIRE_dt   " << FIRE_dt <<std::endl;
//                std::cout << "FIRE_alpha   " << FIRE_alpha  <<std::endl;
//                std::cout << "FIRE_Np+   " << FIRE_N_positive  <<std::endl;
//                std::cout << "FIRE_Np-   " << FIRE_N_negative  <<std::endl;
                std::cout << "\n" << std::endl;
            
            
            if (power >0) {
                 FIRE_N_positive+=1;
                if (FIRE_N_positive> pars.FIRE_N_positive_min){
                    FIRE_dt = fmin(FIRE_dt * pars.FIRE_finc, pars.FIRE_dtmax);
                    FIRE_alpha *= pars.FIRE_falpha;
                }
                
            } else {
                FIRE_dt*=pars.FIRE_fdec;
                FIRE_N_positive =0;
                FIRE_alpha = pars.FIRE_alpha_start;
                vx=0;
                vy=0;
                
            }
                
            
        }
            
    }
