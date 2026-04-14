#include "../include/approximation.hpp"
#include "../include/linear_algebra.hpp"
#include <stdexcept>
#include <cmath>

std::vector<double> polynomialApproximation(const std::vector<double>& x, const std::vector<double>& y, int degree) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Input vectors x and y must have the same size.");
    }
    int n = x.size();
    int m = degree + 1;
    if (n < m) {
        throw std::invalid_argument("Number of points must be at least degree + 1.");
    }

    std::vector<std::vector<double>> A(m, std::vector<double>(m, 0.0));
    std::vector<double> B(m, 0.0);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            double sum = 0;
            for (int k = 0; k < n; ++k) {
                sum += std::pow(x[k], i + j);
            }
            A[i][j] = sum;
        }
        double sum = 0;
        for (int k = 0; k < n; ++k) {
            sum += y[k] * std::pow(x[k], i);
        }
        B[i] = sum;
    }

    auto coeffs = gaussianElimination(A, B);
    if (coeffs) {
        return *coeffs;
    } else {
        throw std::runtime_error("Could not solve the system for approximation coefficients. The system may be ill-conditioned.");
    }
}

double evaluatePolynomial(const std::vector<double>& coeffs, double x) {
    double result = 0.0;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * std::pow(x, i);
    }
    return result;
}