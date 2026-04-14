# Moja Biblioteka Numeryczna

Projekt na przedmiot "Metody Numeryczne". Jest to biblioteka C++ implementująca różne algorytmy numeryczne.

## Funkcjonalność

*   Rozwiązywanie układów równań liniowych (eliminacja Gaussa, dekompozycja LU)
*   Interpolacja (Lagrange, Newton)
*   Aproksymacja
*   Całkowanie numeryczne
*   Rozwiązywanie równań różniczkowych
*   Rozwiązywanie równań nieliniowych

## Jak zbudować i uruchomić (używając VS Code z CMake Tools)

1.  Otwórz folder z projektem w Visual Studio Code.
2.  Upewnij się, że masz zainstalowane rozszerzenia C/C++ i CMake Tools.
3.  VS Code powinien automatycznie skonfigurować projekt. Wybierz kompilator (np. z Visual Studio Build Tools).
4.  Na dolnym pasku stanu:
    *   Wybierz cel do uruchomienia: `[ExamplesApp]` lub `[TestsApp]`.
    *   Kliknij przycisk **Build**, aby skompilować projekt.
    *   Kliknij przycisk **Run (▶️)**, aby uruchomić wybrany program.


## Algebra liniowa (`linear_algebra.hpp`)

### `gaussianElimination(A, b)`
**Opis:** Rozwiązuje układ równań liniowych `Ax = b` metodą eliminacji Gaussa.  
**Argumenty:**
- `A`: macierz współczynników `n×n` (`std::vector<std::vector<double>>`)
- `b`: wektor wyrazów wolnych (`std::vector<double>`)  
**Zwraca:** `std::optional<std::vector<double>>` – wektor rozwiązania `x`, lub `nullopt` jeśli układ jest sprzeczny lub osobliwy.

### `luDecomposition(A)`
**Opis:** Wykonuje dekompozycję LU macierzy `A`, tzn. `A = L·U`.  
**Argumenty:**  
- `A`: macierz `n×n` (`std::vector<std::vector<double>>`)  
**Zwraca:** para macierzy `L` i `U`.

### `forwardSubstitution(L, b)`
**Opis:** Rozwiązuje `Ly = b` dla dolnej trójkątnej macierzy `L`.  
**Zwraca:** `std::vector<double> y`.

### `backwardSubstitution(U, y)`
**Opis:** Rozwiązuje `Ux = y` dla górnej trójkątnej macierzy `U`.  
**Zwraca:** `std::vector<double> x`.


## Interpolacja (`interpolation.hpp`)

### `lagrangeInterpolation(x, y, xp)`
**Opis:** Oblicza wartość interpolowaną w punkcie `xp` metodą Lagrange’a.  
**Argumenty:**  
- `x`, `y`: wektory węzłów interpolacji (`std::vector<double>`)  
- `xp`: punkt, w którym interpolujemy (`double`)  
**Zwraca:** wartość funkcji w `xp`.

### `calculateDividedDifferences(x, y)`
**Opis:** Oblicza współczynniki dzielonych różnic dla interpolacji Newtona.  
**Zwraca:** wektor współczynników (`std::vector<double>`).

### `newtonInterpolation(x, coeffs, xp)`
**Opis:** Oblicza wartość interpolowaną w punkcie `xp` metodą Newtona.  
**Zwraca:** wartość funkcji w `xp`.


## Aproksymacja (`approximation.hpp`)

### `polynomialApproximation(x, y, degree)`
**Opis:** Oblicza współczynniki wielomianu aproksymującego dane metodą najmniejszych kwadratów.  
**Argumenty:**  
- `x`, `y`: dane wejściowe (`std::vector<double>`)  
- `degree`: stopień wielomianu (`int`)  
**Zwraca:** wektor współczynników wielomianu (`std::vector<double>`).

### `evaluatePolynomial(coeffs, x)`
**Opis:** Oblicza wartość wielomianu dla danego `x`.  
**Zwraca:** `double` – wartość funkcji.


## ∫ Całkowanie numeryczne (`integration.hpp`)

### `rectangleMethod(f, a, b, n)`
**Opis:** Całkuje funkcję `f` na przedziale `[a, b]` metodą prostokątów.  
**Zwraca:** przybliżoną wartość całki (`double`).

### `trapezoidalMethod(f, a, b, n)`
**Opis:** Całkuje funkcję `f` na przedziale `[a, b]` metodą trapezów.  
**Zwraca:** `double`.

### `simpsonMethod(f, a, b, n)`
**Opis:** Całkuje funkcję `f` metodą Simpsona (n parzyste).  
**Zwraca:** `double`.

### `compositeGaussLegendre(f, a, b, n_points, subdivisions)`
**Opis:** Stosuje złożoną kwadraturę Gaussa-Legendre’a na przedziale `[a, b]`.  
**Zwraca:** `double`.


## Równania różniczkowe zwyczajne (`differential_equations.hpp`)

Każda funkcja zwraca `std::vector<std::pair<double, double>>`, gdzie pierwszy element to `t`, a drugi to `y(t)`.

### `eulerMethod(f, t0, y0, t_end, h)`
**Opis:** Rozwiązuje ODE metodą Eulera.  
**Zwraca:** wektor punktów rozwiązania.

### `heunMethod(f, t0, y0, t_end, h)`
**Opis:** Rozwiązuje ODE metodą Heuna (prostą predyktor-korektor).  
**Zwraca:** wektor punktów rozwiązania.

### `rk4Method(f, t0, y0, t_end, h)`
**Opis:** Rozwiązuje ODE metodą Rungego-Kutty 4 rzędu.  
**Zwraca:** wektor punktów rozwiązania.


## Równania nieliniowe (`nonlinear_equations.hpp`)

Każda metoda zwraca `std::optional<double>` – wartość pierwiastka lub `nullopt` gdy nie znaleziono.

### `bisection(f, a, b, tol=1e-9, max_iter=1000)`
**Opis:** Znajduje pierwiastek funkcji metodą bisekcji.

### `newtonMethod(f, df, x0, tol=1e-9, max_iter=1000)`
**Opis:** Metoda Newtona z analityczną pochodną.

### `secantMethod(f, x0, x1, tol=1e-9, max_iter=1000)`
**Opis:** Metoda siecznych (bez potrzeby znajomości pochodnej).

### `regulaFalsi(f, a, b, tol=1e-9, max_iter=1000)`
**Opis:** Metoda fałszywej pozycji (Regula Falsi).


## Przykład użycia

```cpp
#include <iostream>
#include "linear_algebra.hpp"

int main() {
    std::vector<std::vector<double>> a = {{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}};
    std::vector<double> b = {8, -11, -3};
    
    auto result = gaussianElimination(a, b);
    if (result) {
        
    }
    return 0;
}

#include "interpolation.hpp"

std::vector<double> x_nodes = {0, 1, 2};
std::vector<double> y_nodes = {0, 1, 4}; 
double interpolated_value = lagrangeInterpolation(x_nodes, y_nodes, 1.5);


#include "approximation.hpp"
// ...
std::vector<double> x = {0, 1, 2};
std::vector<double> y = {1.1, 2.9, 5.2}; 
auto coeffs = polynomialApproximation(x, y, 1);



#include "integration.hpp"
// ...
auto func_to_integrate = [](double x){ return 3 * x * x; };
double integral = simpsonMethod(func_to_integrate, 0.0, 2.0, 100);



#include "differential_equations.hpp"

auto ode_func = [](double t, double y){ return y; };
auto solution_points = rk4Method(ode_func, 0, 1, 1.0, 0.1);
double final_y = solution_points.back().second;



#include "nonlinear_equations.hpp"

auto nonlinear_func = [](double x){ return x * x - 2.0; };
auto root = bisection(nonlinear_func, 1.0, 2.0);
if (root) {
   
}