#pragma once
#include <vector>
#include <array>


struct GlobalData {
    double simulationTime;
    double simulationStepTime;
    double conductivity;
    double alfa;
    double tot;
    double initialTemp;
    double density;
    double specificHeat;
    int nodesNumber;
    int elementsNumber;
    unsigned int npc;

};

struct Node {
    int id;
    double x, y;
    bool BC;
};

struct Jakobian {
    double J[2][2];
    double invJ[2][2];   
    double detJ;
};

struct Element {
    int id;
    int nodeIDs[4]; 
    Jakobian *j;
    std::array<double, 4> P_local;
    std::array<std::array<double, 4>, 4> Hl;
    std::array<std::array<double, 4>, 4> Hbc;
    std::array<std::array<double, 4>, 4> C;
    Element(unsigned int npc) 
    {
        j = new Jakobian[npc * npc];
        for (int i = 0; i < 4; ++i) {
            nodeIDs[i] = 0;
        }
       id = 0;
    }
    Element() {
        j = nullptr;
        for (int i = 0; i < 4; ++i) {
            nodeIDs[i] = 0;
        }
        id = 0;
    }
};

struct Grid {
    int nN; // Liczba węzłów
    int nE; // Liczba elementów
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<int> boundaryConditions; 
    std::vector<std::vector<double>> H_global;
    std::vector<double> P_global;
    std::vector<std::vector<double>> C_global;
};