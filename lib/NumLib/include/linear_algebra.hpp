#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#include <vector>
#include <optional>
#include <utility>

std::optional<std::vector<double>> gaussianElimination(std::vector<std::vector<double>> a, std::vector<double> b);
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> luDecomposition(const std::vector<std::vector<double>>& a);
std::vector<double> forwardSubstitution(const std::vector<std::vector<double>>& L, const std::vector<double>& b);
std::vector<double> backwardSubstitution(const std::vector<std::vector<double>>& U, const std::vector<double>& y);

#endif