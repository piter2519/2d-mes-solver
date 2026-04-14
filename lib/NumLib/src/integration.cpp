#include "../include/integration.hpp"
#include <cmath>
#include <vector>
#include <stdexcept>

double rectangleMethod(std::function<double(double)> f, double a, double b, int n) {
    double h = (b - a) / n;
    double integral = 0.0;
    for (int i = 0; i < n; i++) {
        integral += f(a + i * h);
    }
    return integral * h;
}

double trapezoidalMethod(std::function<double(double)> f, double a, double b, int n) {
    double h = (b - a) / n;
    double integral = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; i++) {
        integral += f(a + i * h);
    }
    return integral * h;
}

double simpsonMethod(std::function<double(double)> f, double a, double b, int n) {
    if (n % 2 != 0) n++; // Simpson's rule requires an even number of intervals
    double h = (b - a) / n;
    double integral = f(a) + f(b);
    for (int i = 1; i < n; i++) {
        integral += (i % 2 == 0 ? 2 : 4) * f(a + i * h);
    }
    return integral * h / 3.0;
}

namespace { // Helper function, not exposed in header
    struct GaussLegendreData {
        int n;
        std::vector<double> nodes;
        std::vector<double> weights;
    };

    GaussLegendreData getGLData(int n) {
        GaussLegendreData data;
        data.n = n;
        if (n == 2) {
            data.nodes = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
            data.weights = {1.0, 1.0};
        } else if (n == 3) {
            data.nodes = {-std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0)};
            data.weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
        } else if (n == 4) {
            data.nodes = {-std::sqrt(3.0/7.0 + 2.0/7.0 * std::sqrt(6.0/5.0)), -std::sqrt(3.0/7.0 - 2.0/7.0 * std::sqrt(6.0/5.0)), std::sqrt(3.0/7.0 - 2.0/7.0 * std::sqrt(6.0/5.0)), std::sqrt(3.0/7.0 + 2.0/7.0 * std::sqrt(6.0/5.0))};
            data.weights = {(18.0 - std::sqrt(30.0))/36.0, (18.0 + std::sqrt(30.0))/36.0, (18.0 + std::sqrt(30.0))/36.0, (18.0 - std::sqrt(30.0))/36.0};
        } else {
            throw std::invalid_argument("Only 2, 3, and 4-point Gauss-Legendre quadrature is supported.");
        }
        return data;
    }
}

double compositeGaussLegendre(std::function<double(double)> f, double a, double b, int n_points, int subdivisions) {
    GaussLegendreData gl_data = getGLData(n_points);
    double total_integral = 0.0;
    double h = (b - a) / subdivisions;

    for (int i = 0; i < subdivisions; ++i) {
        double sub_a = a + i * h;
        double sub_b = a + (i + 1) * h;
        double sub_integral = 0.0;
        for (int j = 0; j < n_points; ++j) {
            double t = gl_data.nodes[j];
            double w = gl_data.weights[j];
            double x = (sub_b - sub_a) / 2.0 * t + (sub_a + sub_b) / 2.0;
            sub_integral += w * f(x);
        }
        total_integral += (sub_b - sub_a) / 2.0 * sub_integral;
    }
    return total_integral;
}