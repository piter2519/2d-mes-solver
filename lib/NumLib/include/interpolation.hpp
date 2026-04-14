#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <vector>

double lagrangeInterpolation(const std::vector<double>& x_nodes, const std::vector<double>& y_nodes, double xp);
double newtonInterpolation(const std::vector<double>& x_nodes, const std::vector<double>& factors, double xp);
std::vector<double> calculateDividedDifferences(const std::vector<double>& x_nodes, const std::vector<double>& y_nodes);

#endif