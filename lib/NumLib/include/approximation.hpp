#ifndef APPROXIMATION_HPP
#define APPROXIMATION_HPP

#include <vector>
#include <functional>

std::vector<double> polynomialApproximation(const std::vector<double>& x, const std::vector<double>& y, int degree);
double evaluatePolynomial(const std::vector<double>& coeffs, double x);

#endif