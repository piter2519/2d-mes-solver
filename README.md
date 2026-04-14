# MESProject

Projekt w C++ (C++17) realizujący symulację niestacjonarnego przewodzenia ciepła metodą elementów skończonych (MES) dla siatki 2D.

Program:
- wczytuje dane z pliku wejściowego w katalogu `data/`,
- buduje macierze globalne (`H`, `C`) oraz wektor (`P`),
- rozwiązuje układ równań w kolejnych krokach czasowych,
- wypisuje na stdout temperaturę minimalną i maksymalną w węzłach,
- zapisuje min/max do pliku wynikowego w `data/my_results/`.

## Wymagania

- Windows 10/11
- Kompilator C++ z obsługą C++17:
  - Visual Studio 2022 (MSVC) **albo** MinGW-w64 (g++)
- CMake >= 3.20

Repo zawiera bibliotekę `NumLib` w `lib/NumLib` (budowana jako część projektu CMake).

## Kompilacja (CMake)

Opcje

### Opcja A: Visual Studio 2022 (MSVC)

```bat
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

### Opcja B: MinGW-w64 (g++)

Jeżeli masz zainstalowane MinGW i jest w `PATH`, możesz użyć generatora MinGW:

```bat
cmake -S . -B build -G "MinGW Makefiles"
cmake --build build -j
```

### Opcja C:

Jakikolwiek inny sposób jaki, który poprawnie kompiluje z użyciem CMake np. wtyczki do Visual Studio Code.

### Gdzie jest plik wykonywalny?

Projekt ma ustawione wyjście binarki na katalog główny repozytorium, więc po zbudowaniu pojawi się:
- `MESProject.exe` (w katalogu głównym)

## Uruchamianie

Uruchamiaj program **z katalogu głównego projektu**, ponieważ pliki wejściowe są wskazywane względnie jako `data/<nazwa_pliku>`.

```bat
.\MESProject.exe
```

### Wybór pliku wejściowego i liczby punktów całkowania

Program **nie przyjmuje argumentów CLI** — plik wejściowy oraz liczba punktów całkowania Gaussa są ustawione na stałe w kodzie w [src/main.cpp](src/main.cpp):
- `std::string filename = "Test1_4_4.txt";`
- `unsigned int npc = 4;`

Aby uruchomić symulację dla innych danych:
1. Zmień `filename` na jeden z plików z katalogu `data/` (np. `Test2_4_4_MixGrid.txt`, `Test3_31_31_kwadrat.txt`).
2. Zmień `npc` (np. 2/3/4 — zależnie od tego co chcesz porównać).
3. Przebuduj projekt.

### Wyniki

W trakcie działania program wypisuje min/max temperatury dla kolejnych kroków czasu.
Dodatkowo zapisuje wynik do pliku:

`data/my_results/<plik_wejściowy>_<npc>pc_results.txt`

Format linii w pliku wynikowym:

`step<TAB>time<TAB>minT<TAB>maxT`

## Format danych wejściowych

Parser wczytuje:
- parametry globalne (czas symulacji, krok, materiał, itp.),
- sekcję `*Node` (węzły),
- sekcję `*Element` (elementy),
- sekcję `*BC` (warunki brzegowe – lista ID węzłów).

Przykładowe pliki wejściowe są w katalogu `data/`.
