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
    double r=0;
    double f=0;
    double nx=0;
    double ny=0;
};




struct PowerlawRepulsion compute_powerlaw_replusion_by_segment(double x, double x1, double x2, double y, double y1,double y2, double sigma, double epsilon, double rcut);

struct PowerlawRepulsion walldiscretePL(double x, double x0, double x1, double y, double y0,double y1, double sigma, double epsilon, double rcut, int numwallnodes);

bool isInsideTriangle(float x, float xo, float x1, float x2, float y, float yo,float y1,float y2);

void shear_special_FIRE(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys);
void shear_special_GD(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys);
void shear_special_stepGD(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys);

