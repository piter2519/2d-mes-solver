#include "../include/linear_algebra.hpp"
#include <cmath>
#include <stdexcept>

std::optional<std::vector<double>> gaussianElimination(std::vector<std::vector<double>> a, std::vector<double> b) {
    const double epsilon = 1e-12;
    size_t n = b.size();
    if (a.size() != n || (n > 0 && a[0].size() != n)) {
        throw std::invalid_argument("Matrix and vector dimensions do not match.");
    }
    for (size_t i = 0; i < n; ++i) {
        size_t maxRow = i;
        for (size_t k = i + 1; k < n; ++k) {
            if (std::abs(a[k][i]) > std::abs(a[maxRow][i])) {
                maxRow = k;
            }
        }
        std::swap(a[i], a[maxRow]);
        std::swap(b[i], b[maxRow]);
        if (std::abs(a[i][i]) < epsilon) return std::nullopt;
        for (size_t k = i + 1; k < n; ++k) {
            double factor = a[k][i] / a[i][i];
            for (size_t j = i; j < n; ++j) {
                a[k][j] -= factor * a[i][j];
            }
            b[k] -= factor * b[i];
        }
    }
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (size_t j = i + 1; j < n; ++j) {
            x[i] -= a[i][j] * x[j];
        }
        x[i] /= a[i][i];
    }
    return x;
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> luDecomposition(const std::vector<std::vector<double>>& a) {
    size_t n = a.size();
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n; j++) {
            double sum = 0.0;
            for (size_t k = 0; k < i; k++) sum += L[i][k] * U[k][j];
            U[i][j] = a[i][j] - sum;
        }
        for (size_t j = i; j < n; j++) {
            if (i == j) {
                if (std::abs(U[i][i]) < 1e-12) throw std::runtime_error("Matrix is singular.");
                L[i][i] = 1.0;
            } else {
                double sum = 0.0;
                for (size_t k = 0; k < i; k++) sum += L[j][k] * U[k][i];
                L[j][i] = (a[j][i] - sum) / U[i][i];
            }
        }
    }
    return {L, U};
}

std::vector<double> forwardSubstitution(const std::vector<std::vector<double>>& L, const std::vector<double>& b) {
    int n = b.size();
    std::vector<double> y(n);
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
    }
    return y;
}

std::vector<double> backwardSubstitution(const std::vector<std::vector<double>>& U, const std::vector<double>& y) {
    int n = y.size();
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
    return x;
}