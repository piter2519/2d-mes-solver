#ifndef NONLINEAR_EQUATIONS_HPP
#define NONLINEAR_EQUATIONS_HPP

#include <functional>
#include <optional>
#include <vector>

std::optional<double> bisection(std::function<double(double)> f, double a, double b, double tol = 1e-9, int max_iter = 100);
std::optional<double> newtonMethod(std::function<double(double)> f, std::function<double(double)> df, double x0, double tol = 1e-9, int max_iter = 100);
std::optional<double> secantMethod(std::function<double(double)> f, double x0, double x1, double tol = 1e-9, int max_iter = 100);
std::optional<double> regulaFalsi(std::function<double(double)> f, double a, double b, double tol = 1e-9, int max_iter = 100);

#endif