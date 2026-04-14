#pragma once
#include "globalData.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

void readDataFromFile(const std::string& filename, Grid& grid, GlobalData& globalData);
void printData(const Grid& grid, const GlobalData& globalData);