#include "../include/nonlinear_equations.hpp"
#include <cmath>

std::optional<double> bisection(std::function<double(double)> f, double a, double b, double tol, int max_iter) {
    if (f(a) * f(b) >= 0.0) return std::nullopt;
    double c = a;
    for (int i = 0; i < max_iter; ++i) {
        c = (a + b) / 2.0;
        if (std::abs(f(c)) < tol || (b - a) / 2.0 < tol) return c;
        if (f(c) * f(a) < 0.0) b = c;
        else a = c;
    }
    return c;
}

std::optional<double> newtonMethod(std::function<double(double)> f, std::function<double(double)> df, double x0, double tol, int max_iter) {
    for (int i = 0; i < max_iter; ++i) {
        double fx = f(x0);
        double dfx = df(x0);
        if (std::abs(dfx) < 1e-12) return std::nullopt; // Pochodna bliska zeru
        double x1 = x0 - fx / dfx;
        if (std::abs(x1 - x0) < tol) return x1;
        x0 = x1;
    }
    return x0; // Zwraca ostatnie przybliżenie po max iteracjach
}

std::optional<double> secantMethod(std::function<double(double)> f, double x0, double x1, double tol, int max_iter) {
    for (int i = 0; i < max_iter; ++i) {
        double fx0 = f(x0);
        double fx1 = f(x1);
        if (std::abs(fx1 - fx0) < 1e-12) return std::nullopt;
        double x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        if (std::abs(x2 - x1) < tol) return x2;
        x0 = x1;
        x1 = x2;
    }
    return x1;
}

std::optional<double> regulaFalsi(std::function<double(double)> f, double a, double b, double tol, int max_iter) {
    if (f(a) * f(b) >= 0) return std::nullopt;
    double x = a;
    for (int i = 0; i < max_iter; ++i) {
        double fa = f(a);
        double fb = f(b);
        if (std::abs(fb - fa) < 1e-12) return std::nullopt;
        x = b - fb * (b - a) / (fb - fa);
        double fx = f(x);
        if (std::abs(fx) < tol) return x;
        if (fa * fx < 0.0) b = x;
        else a = x;
    }
    return x;
}