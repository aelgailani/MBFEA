

#ifndef solvers_hpp
#define solvers_hpp

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
#include "BaseSysData.hpp"
#include "Configuration.hpp"

void gd_solver(const BaseSysData& baseData, const Parameters& pars, long timeStep, int stage, Configuration& mainSys);
void fire_solver(const BaseSysData& baseData, const Parameters& pars, long timeStep, int stage, Configuration& mainSys);

void fire_solver_old(const BaseSysData& baseData, const Parameters& pars, long timeStep, int stage, Configuration& mainSys);

#endif /* solvers_hpp */
