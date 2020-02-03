

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

void GD_solver(const BaseSysData& baseData, const Parameters& pars, int timeStep, int stage, Configuration& mainSys);
void FIRE_solver(const BaseSysData& baseData, const Parameters& pars, int timeStep, int stage, Configuration& mainSys);

#endif /* solvers_hpp */
