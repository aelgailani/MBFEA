#include "Parameters.hpp"
#include "Configuration.hpp"
#include "Parameters.hpp"

void compress(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys);
void stepshear(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys);
void continuousshear(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys);
//void affine_axial_shearing(const BaseSysData& baseData, const Parameters& pars,Configuration& mainSys, double strain);
//void affine_compression(const BaseSysData& baseData, const Parameters& pars, Configuration& mainSys, double strain);
//void hold(const BaseSysData& baseData, Configuration& mainSys, const Parameters& pars);
struct PowerlawRepulsion {
    double fx=0;
    double fy=0;
    double energy=0;
    
};

struct PowerlawRepulsion compute_powerlaw_replusion_by_segment(double x, double x1, double x2, double y, double y1,double y2, double sigma, double epsilon, double rcut);

bool isInsideTriangle(float x, float xo, float x1, float x2, float y, float yo,float y1,float y2);

void shear_special_FIRE(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys);
void shear_special_GD(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys);

