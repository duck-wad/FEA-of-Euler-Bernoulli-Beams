#pragma once

#include <iostream>
#include <vector>
#include <cmath>

std::vector<double> GaussSeidel(const std::vector<std::vector<double>>& A, const std::vector<double>& b, const double errorTol = 0.0001, const int maxIter = 500);
