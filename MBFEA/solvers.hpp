

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

void gd_solver(const BaseSysData& baseData, const Parameters& pars, long& timeStep, long& step,  std::string name, Configuration& mainSys, bool dumpStateData, bool surfaceInteractions);
void fire_solver(const BaseSysData& baseData, const Parameters& pars, long timeStep,  std::string name, Configuration& mainSys);
void fire_solver_old(const BaseSysData& baseData, const Parameters& pars, long timeStep,  std::string name, int stage, Configuration& mainSys);
void fire2_solver(const BaseSysData& baseData, const Parameters& pars, long& timeStep , Configuration& mainSys, long step, bool surfaceInteractions);


