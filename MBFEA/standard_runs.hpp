#include "Parameters.hpp"
#include "Configuration.hpp"
#include "Parameters.hpp"

void compress(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys);
void stepshear(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys);
void continuousshear(const BaseSysData& baseData, const Parameters& pars, long timeStep , Configuration& mainSys);
//void affine_axial_shearing(const BaseSysData& baseData, const Parameters& pars,Configuration& mainSys, double strain);
//void affine_compression(const BaseSysData& baseData, const Parameters& pars, Configuration& mainSys, double strain);
//void hold(const BaseSysData& baseData, Configuration& mainSys, const Parameters& pars);
