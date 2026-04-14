#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <stdexcept>

// Dołączamy wszystkie moduły do testowania
#include "../include/linear_algebra.hpp"
#include "../include/interpolation.hpp"
#include "../include/approximation.hpp"
#include "../include/integration.hpp"
#include "../include/differential_equations.hpp"
#include "../include/nonlinear_equations.hpp"

// --- Funkcje pomocnicze do asercji ---
const double TOL = 1e-9;
void assert_equal(double a, double b, double tol = TOL) { assert(std::abs(a - b) < tol); }
void assert_vectors_equal(const std::vector<double>& v1, const std::vector<double>& v2, double tol = TOL) {
    assert(v1.size() == v2.size());
    for (size_t i = 0; i < v1.size(); ++i) assert_equal(v1[i], v2[i], tol);
}
template<typename Func>
void assert_throws(Func&& func) {
    bool exception_thrown = false;
    try { func(); } catch (const std::exception&) { exception_thrown = true; }
    assert(exception_thrown);
}


// --- Grupy testów dla każdego modułu ---

void run_linear_algebra_tests() {
    std::cout << "Testy: Algebra Liniowa... ";
    // 1. gaussianElimination - poprawny
    std::vector<std::vector<double>> a1 = {{2,1,-1},{-3,-1,2},{-2,1,2}};
    std::vector<double> b1 = {8,-11,-3}, expected_x = {2,3,-1};
    auto result1 = gaussianElimination(a1, b1);
    assert(result1.has_value());
    assert_vectors_equal(*result1, expected_x);
    // 2. gaussianElimination - błędny
    std::vector<std::vector<double>> a2 = {{1,1},{1,1}};
    std::vector<double> b2 = {1,2};
    auto result2 = gaussianElimination(a2, b2);
    assert(!result2.has_value());
    std::cout << "OK\n";
}

void run_interpolation_tests() {
    std::cout << "Testy: Interpolacja... ";
    std::vector<double> xn = {0,1,2}, yn = {0,1,4}; // f(x)=x^2
    // 1. lagrangeInterpolation - poprawny
    assert_equal(lagrangeInterpolation(xn, yn, 1.5), 2.25);
    // 2. lagrangeInterpolation - błędny
    std::vector<double> yn_bad = {0,1};
    assert_throws([&](){ lagrangeInterpolation(xn, yn_bad, 1.5); });
    std::cout << "OK\n";
}

void run_approximation_tests() {
    std::cout << "Testy: Aproksymacja... ";
    std::vector<double> x = {0,1,2}, y = {1,3,5}; // y = 2x+1
    // 1. polynomialApproximation - poprawny
    auto coeffs = polynomialApproximation(x, y, 1);
    assert_equal(coeffs[0], 1.0); assert_equal(coeffs[1], 2.0);
    // 2. polynomialApproximation - błędny
    std::vector<double> x_bad = {0,1};
    assert_throws([&](){ polynomialApproximation(x_bad, y, 2); });
    std::cout << "OK\n";
}

void run_integration_tests() {
    std::cout << "Testy: Calkowanie... ";
    auto f = [](double x) { return 2*x; };
    // 1. simpsonMethod - poprawny
    assert_equal(simpsonMethod(f, 0.0, 3.0, 100), 9.0, 1e-6);
    // 2. simpsonMethod - błędny (początek > koniec)
    assert_equal(simpsonMethod(f, 3.0, 0.0, 100), -9.0, 1e-6);
    std::cout << "OK\n";
}

void run_differential_equations_tests() {
    std::cout << "Testy: Rownania rozniczkowe... ";
    auto f = [](double t, double y) { return y; };
    // 1. rk4Method - poprawny
    auto result1 = rk4Method(f, 0, 1, 1.0, 0.1);
    assert(!result1.empty());
    assert_equal(result1.back().second, std::exp(1.0), 1e-5);
    // 2. rk4Method - błędny (niepoprawny krok h)
    assert_throws([&](){ rk4Method(f, 0, 1, 1.0, 0.0); });
    std::cout << "OK\n";
}

void run_nonlinear_equations_tests() {
    std::cout << "Testy: Rownania nieliniowe... ";
    auto f = [](double x) { return x*x - 4; };
    // 1. bisection - poprawny
    auto r1 = bisection(f, 0, 3);
    assert(r1.has_value()); assert_equal(*r1, 2.0);
    // 2. bisection - błędny
    auto r2 = bisection(f, 2.1, 3.0);
    assert(!r2.has_value());
    std::cout << "OK\n";
}


int main() {
    try {
        run_linear_algebra_tests();
        run_interpolation_tests();
        run_approximation_tests();
        run_integration_tests();
        run_differential_equations_tests();
        run_nonlinear_equations_tests();
        std::cout << "\n--- WSZYSTKIE TESTY ZAKONCZONE POMYSLNIE ---\n";
    } catch (const std::exception& e) {
        std::cerr << "\n\n--- WYSTAPIL KRYTYCZNY BLAD PODCZAS TESTOW: " << e.what() << " ---\n";
        return 1;
    }
    return 0;
}