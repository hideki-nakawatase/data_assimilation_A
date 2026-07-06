#pragma once
#include <Eigen/Dense>

double random_normal_distribution();

Eigen::VectorXd random_fluctuation(double min, double max);

int uniform_random_number(int min, int max);

Eigen::VectorXi random_number_vec(int index);
