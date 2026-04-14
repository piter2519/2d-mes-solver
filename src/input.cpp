#include "input.h"

void printData(const Grid& grid, const GlobalData& globalData) 
{
    std::cout << "--- Dane Globalne ---" << std::endl;
    std::cout << "SimulationTime: " << globalData.simulationTime << std::endl;
    std::cout << "SimulationStepTime: " << globalData.simulationStepTime << std::endl;
    std::cout << "Conductivity: " << globalData.conductivity << std::endl;
    std::cout << "Alfa: " << globalData.alfa << std::endl;
    std::cout << "Tot: " << globalData.tot << std::endl;
    std::cout << "InitialTemp: " << globalData.initialTemp << std::endl;
    std::cout << "Density: " << globalData.density << std::endl;
    std::cout << "SpecificHeat: " << globalData.specificHeat << std::endl;
    std::cout << "Liczba wezlow: " << globalData.nodesNumber << std::endl;
    std::cout << "Liczba elementow: " << globalData.elementsNumber << std::endl;
    std::cout << "\n-----------------------\n" << std::endl;

    std::cout << "--- Wezly (" << grid.nN << ") ---" << std::endl;
    for (const auto& node : grid.nodes) {
        std::cout << "ID: " << node.id << ",\tX: " << node.x << ",\tY: " << node.y << std::endl;
    }
    std::cout << "\n-----------------------\n" << std::endl;

    std::cout << "--- Elementy (" << grid.nE << ") ---" << std::endl;
    for (const auto& element : grid.elements) {
        std::cout << "ID: " << element.id << ",\tWezly: {"
                  << element.nodeIDs[0] << ", "
                  << element.nodeIDs[1] << ", "
                  << element.nodeIDs[2] << ", "
                  << element.nodeIDs[3] << "}" << std::endl;
    }
    std::cout << "\n-----------------------\n" << std::endl;

    std::cout << "--- Warunki Brzegowe (BC) ---" << std::endl;
    for (size_t i = 0; i < grid.boundaryConditions.size(); ++i) {
        std::cout << grid.boundaryConditions[i] << (i == grid.boundaryConditions.size() - 1 ? "" : ", ");
    }
    std::cout << std::endl;
}

void readDataFromFile(const std::string& filename, Grid& grid, GlobalData& globalData) 
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Blad: Nie mozna otworzyc pliku " << filename << std::endl;
        return;
    }

    std::string line, key;
    char comma;

    // Wczytanie danych globalnych
    file >> key >> globalData.simulationTime;
    file >> key >> globalData.simulationStepTime;
    file >> key >> globalData.conductivity;
    file >> key >> globalData.alfa;
    file >> key >> globalData.tot;
    file >> key >> globalData.initialTemp;
    file >> key >> globalData.density;
    file >> key >> globalData.specificHeat;
    file >> key >> key >> globalData.nodesNumber;
    file >> key >> key >> globalData.elementsNumber;

    grid.nN = globalData.nodesNumber;
    grid.nE = globalData.elementsNumber;
    grid.nodes.resize(grid.nN);
    grid.elements.resize(grid.nE);
    grid.P_global.resize(grid.nN, 0.0);

    // Wczytanie danych o węzłach
    while (std::getline(file, line) && line.find("*Node") == std::string::npos);
    for (int i = 0; i < grid.nN; ++i) {
        file >> grid.nodes[i].id >> comma >> grid.nodes[i].x >> comma >> grid.nodes[i].y;
    }

    // Wczytanie danych o elementach
    while (std::getline(file, line) && line.find("*Element") == std::string::npos);
    for (int i = 0; i < grid.nE; ++i) {
        file >> grid.elements[i].id >> comma 
             >> grid.elements[i].nodeIDs[0] >> comma
             >> grid.elements[i].nodeIDs[1] >> comma
             >> grid.elements[i].nodeIDs[2] >> comma
             >> grid.elements[i].nodeIDs[3];
    }
    
    // Wczytanie warunków brzegowych
    while (std::getline(file, line) && line.find("*BC") == std::string::npos);
    
    std::string bc_line;
    std::getline(file, bc_line);
    std::stringstream ss(bc_line);
    int nodeId;
    while(ss >> nodeId) {
        grid.boundaryConditions.push_back(nodeId);
        if (ss.peek() == ',') {
            ss.ignore();
        }
    }

    for(auto& node : grid.nodes) node.BC = false; 
    
    for (int id : grid.boundaryConditions)
        if (id > 0 && id <= grid.nN) grid.nodes[id - 1].BC = true; 
    file.close();
}