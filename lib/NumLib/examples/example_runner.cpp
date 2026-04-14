#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>

// Dołączamy wszystkie moduły z naszej biblioteki
#include "../include/linear_algebra.hpp"
#include "../include/interpolation.hpp"
#include "../include/approximation.hpp"
#include "../include/integration.hpp"
#include "../include/differential_equations.hpp"
#include "../include/nonlinear_equations.hpp"

void run_all_examples() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "===== PLIK DEMONSTRACYJNY: Zastosowania biblioteki NumLib =====\n";

    // --- Algebra Liniowa ---
    std::cout << "\n--- KATEGORIA: Algebra Liniowa ---\n";
    // Przykład 1.1: Rozwiązanie prostego układu
    std::vector<std::vector<double>> a1 = {{4, -2, 1}, {1, 1, 1}, {9, 3, 1}};
    std::vector<double> b1 = {1, 7, 3};
    auto res1 = gaussianElimination(a1, b1);
    std::cout << "Przyklad 1.1: Uklad rownan dla paraboli przechodzacej przez punkty.\n  Wynik: a=" << (*res1)[0] << ", b=" << (*res1)[1] << ", c=" << (*res1)[2] << std::endl;
    // Przykład 1.2: Dekompozycja LU
    auto [L, U] = luDecomposition({{2, 1}, {4, 3}});
    std::cout << "Przyklad 1.2: Dekompozycja LU macierzy 2x2.\n  L[1][0]=" << L[1][0] << ", U[0][1]=" << U[0][1] << std::endl;

    // --- Interpolacja ---
    std::cout << "\n--- KATEGORIA: Interpolacja ---\n";
    std::vector<double> ix = {0, 1, 3}, iy = {1, 3, 13}; // f(x) = x^2 + x + 1
    // Przykład 2.1: Lagrange
    double lag_res = lagrangeInterpolation(ix, iy, 2.0);
    std::cout << "Przyklad 2.1: Estymacja wartosci w punkcie x=2 (Lagrange).\n  Wynik: " << lag_res << " (dokladnie: 7.0)\n";
    // Przykład 2.2: Newton
    auto factors = calculateDividedDifferences(ix, iy);
    double new_res = newtonInterpolation(ix, factors, 2.0);
    std::cout << "Przyklad 2.2: Estymacja wartosci w punkcie x=2 (Newton).\n  Wynik: " << new_res << " (dokladnie: 7.0)\n";

    // --- Aproksymacja ---
    std::cout << "\n--- KATEGORIA: Aproksymacja ---\n";
    std::vector<double> ax = {0, 1, 2, 3, 4}, ay = {1.1, 2.8, 4.2, 5.9, 7.8};
    // Przykład 3.1: Aproksymacja liniowa (stopień 1)
    auto c1 = polynomialApproximation(ax, ay, 1);
    std::cout << "Przyklad 3.1: Aproksymacja liniowa dla danych z szumem.\n  Prosta: y = " << c1[1] << "x + " << c1[0] << std::endl;
    // Przykład 3.2: Aproksymacja kwadratowa (stopień 2)
    std::vector<double> ax2 = {0, 1, 2}, ay2 = {1, 3, 7}; // y = x^2+x+1
    auto c2 = polynomialApproximation(ax2, ay2, 2);
    std::cout << "Przyklad 3.2: Aproksymacja kwadratowa dla danych bez szumu.\n  Wspolczynniki: c0=" << c2[0] << ", c1=" << c2[1] << ", c2=" << c2[2] << std::endl;

    // --- Całkowanie ---
    std::cout << "\n--- KATEGORIA: Calkowanie ---\n";
    auto f_int = [](double x){ return 3*x*x; };
    // Przykład 4.1: Metoda trapezów
    double trap_res = trapezoidalMethod(f_int, 0, 2, 100);
    std::cout << "Przyklad 4.1: Calka z 3*x^2 od 0 do 2 (trapezy).\n  Wynik: " << trap_res << " (dokladnie: 8.0)\n";
    // Przykład 4.2: Metoda Simpsona
    double simp_res = simpsonMethod(f_int, 0, 2, 100);
    std::cout << "Przyklad 4.2: Calka z 3*x^2 od 0 do 2 (Simpson).\n  Wynik: " << simp_res << " (dokladnie: 8.0)\n";
    
    // --- Równania Różniczkowe ---
    std::cout << "\n--- KATEGORIA: Rownania Rozniczkowe ---\n";
    auto f_ode = [](double t, double y){ return -y; }; // y'=-y, y(0)=1 => y(t)=e^(-t)
    // Przykład 5.1: Metoda Eulera (mniej dokładna)
    auto euler_res = eulerMethod(f_ode, 0, 1, 1.0, 0.1);
    std::cout << "Przyklad 5.1: Rozwiazanie y'=-y dla t=1 (Euler).\n  Wynik: " << euler_res.back().second << " (dokladnie: " << std::exp(-1.0) << ")\n";
    // Przykład 5.2: Metoda RK4 (bardziej dokładna)
    auto rk4_res = rk4Method(f_ode, 0, 1, 1.0, 0.1);
    std::cout << "Przyklad 5.2: Rozwiazanie y'=-y dla t=1 (RK4).\n  Wynik: " << rk4_res.back().second << " (dokladnie: " << std::exp(-1.0) << ")\n";

    // --- Równania Nieliniowe ---
    std::cout << "\n--- KATEGORIA: Rownania Nieliniowe ---\n";
    auto f_nonlin = [](double x){ return x*x - 5; };
    auto df_nonlin = [](double x){ return 2*x; };
    // Przykład 6.1: Metoda bisekcji
    auto bis_res = bisection(f_nonlin, 2, 3);
    std::cout << "Przyklad 6.1: Pierwiastek x^2-5 (bisekcja).\n  Wynik: " << *bis_res << " (dokladnie: " << std::sqrt(5.0) << ")\n";
    // Przykład 6.2: Metoda Newtona
    auto newt_res = newtonMethod(f_nonlin, df_nonlin, 2.0);
    std::cout << "Przyklad 6.2: Pierwiastek x^2-5 (Newton).\n  Wynik: " << *newt_res << " (dokladnie: " << std::sqrt(5.0) << ")\n";
}

int main() {
    try {
        run_all_examples();
        std::cout << "\n\n--- Wszystkie przyklady zostaly wykonane pomyslnie. ---\n";
    } catch (const std::exception& e) {
        std::cerr << "\n\n--- WYSTAPIL KRYTYCZNY BLAD: " << e.what() << " ---\n";
        return 1;
    }
    return 0;
}