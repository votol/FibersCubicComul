#include "ParametersInit.h"
#include "schema.h"
#include "ParameterDefines.h"
#include <iostream>
ParametersHolder::ParametersHolder(YAML::Node&& in)
{
    parameters.resize(5);
    t = 0.0;
    Nfibs = in[FibersCubicComulSchema::PARAMETER_Nfibs].as<double>();
    Nfibs_calc = in[FibersCubicComulSchema::PARAMETER_Nfibs_calc].as<double>();
    L = in[FibersCubicComulSchema::PARAMETER_L].as<double>();
    gamma = in[FibersCubicComulSchema::PARAMETER_gamma].as<double>();
}
