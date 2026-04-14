#pragma once
#include <array>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <functional>
#include "globalData.h"

class GaussLegendreData  {
    public:
    std::array<double, 4> nodes;
    std::array<double, 4> weights;
    GaussLegendreData(int n);
    GaussLegendreData() {}
};

double gaussLegendreIntegrate(std::function<double(double)> f, double a, double b, int n);
double gaussLegendreIntegrate(std::function<double(double, double)> f, double a, double b, int n);
GaussLegendreData getGLData(int n);

double f1(double x);
double f2(double x, double y);