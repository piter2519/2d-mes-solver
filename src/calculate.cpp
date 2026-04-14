#include "calculate.h"

#define DEBUGMODE false

void computeElementH(Element& el, const Grid& grid, const GlobalData& globalData, const ElemUniv& elemUniv, const GaussLegendreData& gauss)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            el.Hl[i][j] = 0.0;

    const unsigned int npc = globalData.npc;

    if(DEBUGMODE)
        std::cout << "Obliczenia macierzy Jakobiego dla elementu nr " << el.id << "\n";

    for (unsigned int pc = 0; pc < npc * npc; ++pc)
    {
        Jakobian J = calculateJakobian(elemUniv, grid, el, pc);

        if(DEBUGMODE) 
        {
            std::cout << std::setprecision(6);
            std::cout << "Macierz Jakobiego dla " << pc + 1 << " punktu całkowania\n";
            std::cout << std::setw(15) << J.J[0][0] << "  " << std::setw(15) << J.J[0][1] << "\n";
            std::cout << std::setw(15) << J.J[1][0] << "  " << std::setw(15) << J.J[1][1] << "\n";
            std::cout << "detJ = " << std::setprecision(9) << J.detJ << "\n\n";
            std::cout << std::setprecision(2);
        }
        std::array<double, 4> dNdx, dNdy;
        //print elemuniv
        if(DEBUGMODE)
        {
            for(int i = 0; i < 4; ++i) 
            {
                std::cout << "dNdKsi[" << i << "] = " << std::setw(6) << elemUniv.dNdKsi[pc][i] 
                        << ", dNdEta[" << i << "] = " << std::setw(6) << elemUniv.dNdEta[pc][i] << "\n";
            }
        }


        dNdx[0] = J.invJ[0][0] * elemUniv.dNdKsi[pc][0] + J.invJ[1][0] * elemUniv.dNdEta[pc][0];
        dNdx[1] = J.invJ[0][0] * elemUniv.dNdKsi[pc][1] + J.invJ[1][0] * elemUniv.dNdEta[pc][1];
        dNdx[2] = J.invJ[0][0] * elemUniv.dNdKsi[pc][2] + J.invJ[1][0] * elemUniv.dNdEta[pc][2];
        dNdx[3] = J.invJ[0][0] * elemUniv.dNdKsi[pc][3] + J.invJ[1][0] * elemUniv.dNdEta[pc][3];

        dNdy[0] = J.invJ[0][1] * elemUniv.dNdKsi[pc][0] + J.invJ[1][1] * elemUniv.dNdEta[pc][0];
        dNdy[1] = J.invJ[0][1] * elemUniv.dNdKsi[pc][1] + J.invJ[1][1] * elemUniv.dNdEta[pc][1];
        dNdy[2] = J.invJ[0][1] * elemUniv.dNdKsi[pc][2] + J.invJ[1][1] * elemUniv.dNdEta[pc][2];
        dNdy[3] = J.invJ[0][1] * elemUniv.dNdKsi[pc][3] + J.invJ[1][1] * elemUniv.dNdEta[pc][3];

        auto Hl_local = calculateMatrix_Hl(dNdx, dNdy, globalData.conductivity, J.detJ);

        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                el.Hl[i][j] += Hl_local[i][j]
                             * gauss.weights[pc / npc]
                             * gauss.weights[pc % npc];
            }
        }

        el.j[pc] = J;
    }

    if(DEBUGMODE)
    {
        std::cout << "Macierz H dla elementu nr " << el.id << "\n";
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                std::cout << std::setw(15) << el.Hl[i][j] << "\t";
            }
            std::cout << '\n';
        }
        std::cout << "\n------------------------------------\n";
    }
}

void assembleElementToGlobal(const Element& el, Grid& grid)
{
    for (int i = 0; i < 4; ++i)
    {
        const int row = el.nodeIDs[i] - 1;
        for (int j = 0; j < 4; ++j)
        {
            const int col = el.nodeIDs[j] - 1;
            grid.H_global[row][col] += el.Hl[i][j] + el.Hbc[i][j];
        }
    }
    // agregacja do P_global
    for (int i = 0; i < 4; i++)
        grid.P_global[el.nodeIDs[i] - 1] += el.P_local[i];

    // agregacja do C_global
    for (int i = 0; i < 4; i++)
    {
        const int row = el.nodeIDs[i] - 1;
        for (int j = 0; j < 4; j++)
        {
            const int col = el.nodeIDs[j] - 1;
            grid.C_global[row][col] += el.C[i][j];
        }
    }
}
void processAllElements(Grid& grid, GlobalData& globalData, const ElemUniv& elemUniv)
{
    GaussLegendreData gauss(globalData.npc);

    for (auto& el : grid.elements)
    {
        computeElementH(el, grid, globalData, elemUniv, gauss);

        // Hbc + P_local + agregacja P_global
        calculateHbc(grid, globalData, elemUniv, el);

        calculateC(grid, globalData, elemUniv, el);

        // Agregacja Hl + Hbc do macierzy globalnej
        assembleElementToGlobal(el, grid);
    }
}
void printGlobalResults(const Grid& grid, const GlobalData& globalData)
{
    std::cout << "Macierz H globalna:\n";
    for (int i = 0; i < globalData.nodesNumber; ++i)
    {
        for (int j = 0; j < globalData.nodesNumber; ++j)
        {
            std::cout << grid.H_global[i][j] << "\t";
        }
        std::cout << '\n';
    }
    std::cout << "\n------------------------------------\n";

    std::cout << "Wektor P globalny:\n";
    for (int i = 0; i < globalData.nodesNumber; ++i)
    {
        std::cout << grid.P_global[i] << "\t";
    }
    std::cout << '\n';

    std::cout << "\n------------------------------------\n";
    std::cout << "Macierz C globalna:\n";
    for (int i = 0; i < globalData.nodesNumber; ++i)
    {
        for (int j = 0; j < globalData.nodesNumber; ++j)
        {
            std::cout << grid.C_global[i][j] << "\t";
        }
        std::cout << '\n';
    }
}


void allocateJakobians(Grid& grid, unsigned int npc)
{
    const unsigned int totalPC = npc * npc;
    for (auto& el : grid.elements)
    {
        el.j = new Jakobian[totalPC];
    }
}

void initializeGlobalH(Grid& grid)
{
    grid.H_global.assign(grid.nN, std::vector<double>(grid.nN, 0.0));
}

void initializeGlobalC(Grid& grid)
{
    grid.C_global.assign(grid.nN, std::vector<double>(grid.nN, 0.0));
}

void initializeSimulation(const std::string& filename, unsigned int npc, Grid& grid, GlobalData& globalData, ElemUniv& elemUniv)
{
{
    readDataFromFile(filename, grid, globalData);
    globalData.npc = npc;
    allocateJakobians(grid, npc);
    initializeGlobalH(grid);
    initializeGlobalC(grid);
}
}

void calculateC(Grid& grid, const GlobalData& globalData, const ElemUniv& eUniv, Element& el) 
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            el.C[i][j] = 0.0;

    GaussLegendreData gauss(globalData.npc);

    
    double rho = globalData.density;
    double c = globalData.specificHeat;
    for (unsigned int pc = 0; pc < globalData.npc * globalData.npc; ++pc)
    {
        double dV = el.j[pc].detJ;
        double w1 = gauss.weights[pc / globalData.npc];
        double w2 = gauss.weights[pc % globalData.npc];
        double commonC = rho * c * dV * w1 * w2;
        double* N = eUniv.N[pc];
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                el.C[i][j] += commonC * N[i] * N[j];
    }


    //print C
    if(DEBUGMODE)
    {
        std::cout << "Macierz C dla elementu nr " << el.id << ":\n";
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                std::cout << std::setw(6) << el.C[i][j] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << "\n------------------------------------\n";
    }
}

void calculateHbc(Grid& grid, const GlobalData& globalData, const ElemUniv& eUniv, Element& el) 
{

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            el.Hbc[i][j] = 0.0;

    GaussLegendreData gauss(globalData.npc);

    for (int i = 0; i < 4; ++i)
        el.P_local[i] = 0.0;

    for (int e = 0; e < 4; e++)
    {
        
        int globA = el.nodeIDs[e];
        int globB = el.nodeIDs[(e+1) % 4];

        bool isA_BC = grid.nodes[globA - 1].BC;
        bool isB_BC = grid.nodes[globB - 1].BC;

        if (!(isA_BC && isB_BC))
            continue;   // ten bok nie ma konwekcji

        const Node& nodeA = grid.nodes[globA - 1];
        const Node& nodeB = grid.nodes[globB - 1];

        double dx = nodeB.x - nodeA.x;
        double dy = nodeB.y - nodeA.y;
        double l = std::sqrt(dx * dx + dy * dy);
        double detJ = l / 2.0;


        for (unsigned int pc = 0; pc < globalData.npc; ++pc)
        {
            double w = gauss.weights[pc];
            double* N = eUniv.s[e].N[pc];
            double commonH = globalData.alfa * detJ * w;
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    el.Hbc[i][j] += commonH * N[i] * N[j];
        }

        // P_local
        for (unsigned int pc = 0; pc < globalData.npc; ++pc)
        {
            double w = gauss.weights[pc];
            double* N = eUniv.s[e].N[pc];
            double commonP = globalData.alfa * globalData.tot * detJ * w;
            for (int i = 0; i < 4; i++)
                el.P_local[i] += commonP * N[i];
        }
    }
}

std::array<std::array<double, 4>, 4> calculateMatrix_Hl(std::array<double, 4> dNdx, std::array<double, 4> dNdy, double k, double dV)
{
    std::array<std::array<double, 4>, 4> Hl;

    for (int i = 0; i < 4; ++i) 
    {
        for (int j = 0; j < 4; ++j) 
        {
            Hl[i][j] = (dNdx[i] * dNdx[j] + dNdy[i] * dNdy[j]) * k * dV;
        }
    }
    return Hl;
}

Jakobian calculateJakobian(const ElemUniv& elemUniv, const Grid& grid, const Element& el, unsigned int pcIndex) 
{
    Jakobian J{};
    J.J[0][0] = J.J[0][1] = J.J[1][0] = J.J[1][1] = 0.0;

    for (int i = 0; i < 4; ++i) {
        double x = grid.nodes[el.nodeIDs[i] - 1].x;
        double y = grid.nodes[el.nodeIDs[i] - 1].y;

        J.J[0][0] += elemUniv.dNdKsi[pcIndex][i] * x;
        J.J[0][1] += elemUniv.dNdEta[pcIndex][i] * x;
        J.J[1][0] += elemUniv.dNdKsi[pcIndex][i] * y;
        J.J[1][1] += elemUniv.dNdEta[pcIndex][i] * y;   
    }

    J.detJ = J.J[0][0] * J.J[1][1] - J.J[0][1] * J.J[1][0];

    double invDet = 1.0 / J.detJ;
    J.invJ[0][0] =  J.J[1][1] * invDet;
    J.invJ[0][1] = -J.J[0][1] * invDet;
    J.invJ[1][0] = -J.J[1][0] * invDet;
    J.invJ[1][1] =  J.J[0][0] * invDet;

    // if(DEBUGMODE)
    // {
    //     std::cout << "\nElement " << el.id << ", PC " << pcIndex + 1 << "\n";
    //     std::cout << "Macierz J:\n";
    //     std::cout << "[" << J.J[0][0] << ", " << J.J[0][1] << "]\n";
    //     std::cout << "[" << J.J[1][0] << ", " << J.J[1][1] << "]\n";
    //     std::cout << "Wyznacznik detJ = " << J.detJ << "\n";
    //     std::cout << "Macierz odwrotna J^-1:\n";
    //     std::cout << "[" << J.invJ[0][0] << ", " << J.invJ[0][1] << "]\n";
    //     std::cout << "[" << J.invJ[1][0] << ", " << J.invJ[1][1] << "]\n";
    //     std::cout << "------------------------------------\n";
    // }
    return J;
}
