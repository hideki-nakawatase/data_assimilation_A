#pragma once
#include <vector>

std::vector<double> lorenz96_rk4(int N, int steps, double dt, double F, std::vector<double> &X);
