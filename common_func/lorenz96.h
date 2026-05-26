#pragma once
#include <Eigen/Dense>

Eigen::VectorXd lorenz96(Eigen::VectorXd &X, int N, double F);
