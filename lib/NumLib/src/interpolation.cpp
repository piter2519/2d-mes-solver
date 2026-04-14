#include "../include/interpolation.hpp"
#include <stdexcept>

double lagrangeInterpolation(const std::vector<double>& x, const std::vector<double>& y, double xp) {
    if (x.size() != y.size() || x.empty()) {
        throw std::invalid_argument("Node vectors must have the same, non-zero size.");
    }
    double yp = 0;
    size_t n = x.size();
    for (size_t i = 0; i < n; i++) {
        double p = y[i];
        for (size_t j = 0; j < n; j++) {
            if (i != j) {
                if (std::abs(x[i] - x[j]) < 1e-12) {
                    throw std::runtime_error("Duplicate x nodes detected, interpolation failed.");
                }
                p *= (xp - x[j]) / (x[i] - x[j]);
            }
        }
        yp += p;
    }
    return yp;
}

std::vector<double> calculateDividedDifferences(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.empty()) {
        throw std::invalid_argument("Node vectors must have the same, non-zero size.");
    }
    size_t n = x.size();
    std::vector<std::vector<double>> table(n, std::vector<double>(n));
    for (size_t i = 0; i < n; ++i) table[i][0] = y[i];
    for (size_t j = 1; j < n; ++j) {
        for (size_t i = j; i < n; ++i) {
            if (std::abs(x[i] - x[i - j]) < 1e-12) {
                throw std::runtime_error("Duplicate x nodes detected, cannot calculate divided differences.");
            }
            table[i][j] = (table[i][j - 1] - table[i - 1][j - 1]) / (x[i] - x[i - j]);
        }
    }
    std::vector<double> factors(n);
    for (size_t i = 0; i < n; ++i) factors[i] = table[i][i];
    return factors;
}

double newtonInterpolation(const std::vector<double>& x, const std::vector<double>& factors, double xp) {
    double result = 0;
    double term = 1;
    for (size_t i = 0; i < x.size(); ++i) {
        result += factors[i] * term;
        term *= (xp - x[i]);
    }
    return result;
}