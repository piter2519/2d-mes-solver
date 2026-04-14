#ifndef INTEGRATION_HPP
#define INTEGRATION_HPP

#include <functional>

double rectangleMethod(std::function<double(double)> f, double a, double b, int n);
double trapezoidalMethod(std::function<double(double)> f, double a, double b, int n);
double simpsonMethod(std::function<double(double)> f, double a, double b, int n);
double compositeGaussLegendre(std::function<double(double)> f, double a, double b, int n_points, int subdivisions);

#endif