#pragma once
#include <Eigen/Dense>

Eigen::VectorXd lorenz96_rk4(int steps, Eigen::VectorXd &X);
