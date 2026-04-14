#ifndef DIFFERENTIAL_EQUATIONS_HPP
#define DIFFERENTIAL_EQUATIONS_HPP

#include <vector>
#include <utility>
#include <functional>

using ODEResult = std::vector<std::pair<double, double>>;

ODEResult eulerMethod(std::function<double(double, double)> f, double t0, double y0, double t_end, double h);
ODEResult heunMethod(std::function<double(double, double)> f, double t0, double y0, double t_end, double h);
ODEResult rk4Method(std::function<double(double, double)> f, double t0, double y0, double t_end, double h);

#endif