#include "../include/differential_equations.hpp"
#include <stdexcept>
#include <cmath>

ODEResult eulerMethod(std::function<double(double, double)> f, double t0, double y0, double t_end, double h) {
    if (h <= 0) {
        throw std::invalid_argument("Step h must be positive.");
    }
    ODEResult results;
    double t = t0, y = y0;
    results.emplace_back(t, y);
    
    int num_steps = static_cast<int>((t_end - t0) / h);
    for (int i = 0; i < num_steps; ++i) {
        y += h * f(t, y);
        t += h;
        results.emplace_back(t, y);
    }
    return results;
}

ODEResult heunMethod(std::function<double(double, double)> f, double t0, double y0, double t_end, double h) {
    if (h <= 0) {
        throw std::invalid_argument("Step h must be positive.");
    }
    ODEResult results;
    double t = t0, y = y0;
    results.emplace_back(t, y);
    
    int num_steps = static_cast<int>((t_end - t0) / h);
    for (int i = 0; i < num_steps; ++i) {
        double k1 = f(t, y);
        double k2 = f(t + h, y + h * k1);
        y += h * 0.5 * (k1 + k2);
        t += h;
        results.emplace_back(t, y);
    }
    return results;
}

ODEResult rk4Method(std::function<double(double, double)> f, double t0, double y0, double t_end, double h) {
    if (h <= 0) {
        throw std::invalid_argument("Step h must be positive.");
    }
    ODEResult results;
    double t = t0, y = y0;
    results.emplace_back(t, y);
    
    int num_steps = static_cast<int>((t_end - t0) / h);
    for (int i = 0; i < num_steps; ++i) {
        double k1 = f(t, y);
        double k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
        double k3 = f(t + 0.5 * h, y + 0.5 * h * k2);
        double k4 = f(t + h, y + h * k3);
        y += (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        t += h;
        results.emplace_back(t, y);
    }
    return results;
}