#include "gauss.h"

GaussLegendreData::GaussLegendreData(int n)
{
    switch (n) 
    {
        case 2:
            this->nodes[0] = -1.0 / std::sqrt(3);
            this->nodes[1] =  1.0 / std::sqrt(3);
            this->weights[0] = 1.0;
            this->weights[1] = 1.0;
            break;
        case 3:
            this->nodes[0] = -std::sqrt(3.0 / 5.0);
            this->nodes[1] = 0.0;
            this->nodes[2] = std::sqrt(3.0 / 5.0);
            this->weights[0] = 5.0 / 9.0;
            this->weights[1] = 8.0 / 9.0;
            this->weights[2] = 5.0 / 9.0;
            break;
        case 4:
            this->nodes[0] = -std::sqrt((3 + 2 * std::sqrt(6.0 / 5.0)) / 7);
            this->nodes[1] = -std::sqrt((3 - 2 * std::sqrt(6.0 / 5.0)) / 7);
            this->nodes[2] =  std::sqrt((3 - 2 * std::sqrt(6.0 / 5.0)) / 7);
            this->nodes[3] =  std::sqrt((3 + 2 * std::sqrt(6.0 / 5.0)) / 7);
            this->weights[0] = (18 - std::sqrt(30)) / 36;
            this->weights[1] = (18 + std::sqrt(30)) / 36;
            this->weights[2] = (18 + std::sqrt(30)) / 36;
            this->weights[3] = (18 - std::sqrt(30)) / 36;
            break;
        default:
            std::cerr << "Obsługiwane są tylko kwadratury 2, 3 i 4 punktowe." << std::endl;
            std::exit(1);
    }
}

double gaussLegendreIntegrate(std::function<double(double)> f, double a, double b, int n) 
{
    GaussLegendreData data(n);
    double result = 0.0;
    for (int i = 0; i < n; ++i) 
    {
        double t = data.nodes[i];
        double x = (b - a) / 2 * t + (a + b) / 2; 
        result += data.weights[i] * f(x);
    }
    result *= (b - a) / 2;
    return result;
}

double gaussLegendreIntegrate(std::function<double(double, double)> f, double a, double b, int n) 
{   
    GaussLegendreData data(n);
    double result = 0.0;
    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < n; ++j) 
        {
            double tx = data.nodes[i];
            double ty = data.nodes[j];
            double x = (b - a) / 2 * tx + (a + b) / 2;
            double y = (b - a) / 2 * ty + (a + b) / 2; 
            result += data.weights[i] * data.weights[j] * f(x, y);
        }
    }
    result *= (b - a) / 2;
    result *= (b - a) / 2;
    return result;
}

GaussLegendreData getGLData(int n) 
{
    GaussLegendreData data;
    switch (n) 
    {
        case 2:
            data.nodes[0] = -1.0 / std::sqrt(3);
            data.nodes[1] =  1.0 / std::sqrt(3);
            data.weights[0] = 1.0;
            data.weights[1] = 1.0;
            break;
        case 3:
            data.nodes[0] = -std::sqrt(3.0 / 5.0);
            data.nodes[1] = 0.0;
            data.nodes[2] = std::sqrt(3.0 / 5.0);
            data.weights[0] = 5.0 / 9.0;
            data.weights[1] = 8.0 / 9.0;
            data.weights[2] = 5.0 / 9.0;
            break;
        case 4:
            data.nodes[0] = -std::sqrt((3 + 2 * std::sqrt(6.0 / 5.0)) / 7);
            data.nodes[1] = -std::sqrt((3 - 2 * std::sqrt(6.0 / 5.0)) / 7);
            data.nodes[2] =  std::sqrt((3 - 2 * std::sqrt(6.0 / 5.0)) / 7);
            data.nodes[3] =  std::sqrt((3 + 2 * std::sqrt(6.0 / 5.0)) / 7);
            data.weights[0] = (18 - std::sqrt(30)) / 36;
            data.weights[1] = (18 + std::sqrt(30)) / 36;
            data.weights[2] = (18 + std::sqrt(30)) / 36;
            data.weights[3] = (18 - std::sqrt(30)) / 36;
            break;
        default:
            std::cerr << "Obsługiwane są tylko kwadratury 2, 3 i 4 punktowe." << std::endl;
            std::exit(1);
    }
    return data;
}

double f1(double x) 
{
    return 5 * x * x + 3 * x + 6;
}

double f2(double x, double y) 
{
    return 5 * x * x * y * y + 3 * x * y + 6;
}
