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


struct closestPointOnHermitianCurve {
    double g=0;
    double x=0;
    double y=0;
    double nx=0;
    double ny=0;
    double gap=0;
    
};

struct point {
    double x=0;
    double y=0;
};


struct closestPointOnHermitianCurve find_closest_point_on_Hermitian_interpolation(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4, double alpha);

struct point HermitianInterpolation(double x2, double x3, double x4, double y2, double y3, double y4, double alpha, double g);


void testFire2(const Parameters& pars);
void testFire(const Parameters& pars);
