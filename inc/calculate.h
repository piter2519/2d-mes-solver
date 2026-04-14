#pragma once
#include <iomanip>
#include "globalData.h"
#include "gauss.h"
#include "elemUniv.h"
#include "input.h"

Jakobian calculateJakobian(const ElemUniv& elemUniv, const Grid& grid, const Element& el, unsigned int pcIndex);
std::array<std::array<double, 4>, 4> calculateMatrix_Hl(std::array<double, 4> dNdx, std::array<double, 4> dNdy, double k, double dV);
void calculateHbc(Grid& grid, const GlobalData& globalData, const ElemUniv& eUniv, Element& el);
void initializeSimulation(const std::string& filename, unsigned int npc, Grid& grid, GlobalData& globalData, ElemUniv& elemUniv);
void allocateJakobians(Grid& grid, unsigned int npc);
void computeElementH(Element& el, const Grid& grid, const GlobalData& globalData, const ElemUniv& elemUniv, const GaussLegendreData& gauss);
void assembleElementToGlobal(const Element& el, Grid& grid);
void processAllElements(Grid& grid, GlobalData& globalData, const ElemUniv& elemUniv);
void printGlobalResults(const Grid& grid, const GlobalData& globalData);
void calculateC(Grid& grid, const GlobalData& globalData, const ElemUniv& eUniv, Element& el);