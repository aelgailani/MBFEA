#include "Configuration.hpp"

// integrator number 0
void semi_implicit_Euler(Configuration& mainSys, double dt, bool FIRE=false, double power=0, double scale1=0, double scale2=0){

    mainSys.velocityX += mainSys.forceX* dt;
    mainSys.velocityY += mainSys.forceY* dt;
    
    if (FIRE){
        if (power>0.0){
            mainSys.velocityX = scale1 * mainSys.velocityX+ scale2 * mainSys.forceX;
            mainSys.velocityY = scale1 * mainSys.velocityY+ scale2 * mainSys.forceY;
        }
    }
    mainSys.curPosX += mainSys.velocityX* dt;
    mainSys.curPosY += mainSys.velocityY* dt;
}

// integrator number 1
void leapfrog(Configuration& mainSys, double dt, bool FIRE=false, double power=0, double scale1=0, double scale2=0){

    mainSys.velocityX += mainSys.forceX* dt;
    mainSys.velocityY += mainSys.forceY* dt;
    
    if (FIRE){
        if (power>0.0){
            mainSys.velocityX = scale1 * mainSys.velocityX+ scale2 * mainSys.forceX;
            mainSys.velocityY = scale1 * mainSys.velocityY+ scale2 * mainSys.forceY;
        }
    }
    mainSys.curPosX += mainSys.velocityX* dt;
    mainSys.curPosY += mainSys.velocityY* dt;
}


// integrator number 2
void explicit_Euler(Configuration& mainSys, double dt, bool FIRE=false, double power=0, double scale1=0, double scale2=0){
    
    if (FIRE){
        if (power>0.0){
            mainSys.velocityX = scale1 * mainSys.velocityX+ scale2 * mainSys.forceX;
            mainSys.velocityY = scale1 * mainSys.velocityY+ scale2 * mainSys.forceY;
        }
        mainSys.curPosX += mainSys.velocityX* dt;
        mainSys.curPosY += mainSys.velocityY* dt;
        mainSys.velocityX += mainSys.forceX* dt;
        mainSys.velocityY += mainSys.forceY* dt;
    }else{

        mainSys.curPosX += mainSys.forceX* dt;
        mainSys.curPosX += mainSys.forceY* dt;
    }
   
    
}


// integrator number 3
void velocity_Verlet(Configuration& mainSys, double dt, bool FIRE=false, double power=0, double scale1=0, double scale2=0){
    
    mainSys.velocityX += mainSys.forceX* dt;
    mainSys.velocityY += mainSys.forceY* dt;
    
    if (FIRE){
        if (power>0.0){
            mainSys.velocityX = scale1 * mainSys.velocityX+ scale2 * mainSys.forceX;
            mainSys.velocityY = scale1 * mainSys.velocityY+ scale2 * mainSys.forceY;
        }
    }
    mainSys.curPosX += mainSys.velocityX* dt;
    mainSys.curPosY += mainSys.velocityY* dt;
    
    mainSys.velocityX += mainSys.forceX* dt;
    mainSys.velocityY += mainSys.forceY* dt;
}
