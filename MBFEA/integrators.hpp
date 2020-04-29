
#include "Configuration.hpp"

void semi_implicit_Euler(Configuration& mainSys, double dt, bool FIRE=false, double power=0, double scale1=0, double scale2=0);
void leapfrog(Configuration& mainSys, double dt, bool FIRE=false, double power=0, double scale1=0, double scale2=0);
void explicit_Euler(Configuration& mainSys, double dt, bool FIRE=false, double power=0, double scale1=0, double scale2=0);
void velocity_Verlet(Configuration& mainSys, double dt, bool FIRE=false, double power=0, double scale1=0, double scale2=0);
