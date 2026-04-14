#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <array>
#include <functional>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "../inc/calculate.h"
#include "linear_algebra.hpp"

void printMinMaxTemp(int step, double time, const std::vector<double>& T) 
{
    double minT = T[0];
    double maxT = T[0];
    for (double val : T) {
        if (val < minT) minT = val;
        if (val > maxT) maxT = val;
    }
    std::cout << "Time: " << std::setw(6) << time << " s (Step " << step << ") | "
              << "Min: " << minT << " | Max: " << maxT << std::endl;
}

void saveMinMaxTempToFile(std::ofstream& file, int step, double time, const std::vector<double>& T) 
{
    double minT = T[0];
    double maxT = T[0];
    for (double val : T) {
        if (val < minT) minT = val;
        if (val > maxT) maxT = val;
    }
    file << std::setprecision(16);
    file << step << "\t" << time << "\t" << minT << "\t" << maxT << "\n";
}

int main() 
{
    GlobalData globalData;
    Grid grid;
    // std::string filename = "Test3_31_31_kwadrat.txt";
    // std::string filename = "Test2_4_4_MixGrid.txt";
    std::string filename = "Test1_4_4.txt";
    unsigned int npc = 4;
    globalData.npc = npc;
    ElemUniv elemUniv(npc);

    std::string pathFilename = "data/" + filename;

    std::cout << std::fixed << std::setprecision(2);
    initializeSimulation(pathFilename, npc, grid, globalData, elemUniv);
    processAllElements(grid, globalData, elemUniv);


    double dt = globalData.simulationStepTime;
    double totalTime = globalData.simulationTime;
    
    std::vector<double> t0(globalData.nodesNumber, globalData.initialTemp);
    
    std::vector<std::vector<double>> H_eff = grid.H_global;
    for (int i = 0; i < grid.nN; ++i) {
        for (int j = 0; j < grid.nN; ++j) {
            H_eff[i][j] += grid.C_global[i][j] / dt;
        }
    }

    std::cout << "Rozpoczynam symulacje niestacjonarna..." << std::endl;
    std::cout << "Czas calkowity: " << totalTime << ", Krok: " << dt << std::endl;

    std::string resultFilename = "data/my_results/" + filename + '_' + std::to_string(npc) + "pc_results.txt";
    std::ofstream resultFile(resultFilename, std::ios::out);
    if (!resultFile.is_open()) 
        std::cout << "Nie mozna otworzyc pliku do zapisu wynikow: " << resultFilename << std::endl;

    int step = 0;
    for (double currentTime = dt; currentTime <= totalTime; currentTime += dt) 
    {
        step++;
        std::vector<double> P_eff = grid.P_global;
        
        for (int i = 0; i < grid.nN; ++i) {
            for (int j = 0; j < grid.nN; ++j) {
                P_eff[i] += (grid.C_global[i][j] / dt) * t0[j];
            }
        }
        
        // std::vector<std::vector<double>> A = H_eff;
        // std::vector<double> B = P_eff;
        
        auto result = gaussianElimination(H_eff, P_eff);
        
        if (!result) {
            std::cerr << "Blad rozwiazywania ukladu rownan w kroku " << step << std::endl;
            break;
        }

        t0 = *result;

        printMinMaxTemp(step, currentTime, t0);

        if (resultFile.is_open()) 
            saveMinMaxTempToFile(resultFile, step, currentTime, t0);

        // // temperatura w kazdym z wezlow
        // for(size_t i = 0; i < t0.size(); ++i) 
        // {
        //     std::cout << "T0[" << i + 1 << "] = " << t0[i] << "\n";
        // }
    }
    // printGlobalResults(grid, globalData);
    
    return 0;
}