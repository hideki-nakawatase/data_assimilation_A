#pragma once
#include <Eigen/Dense>

Eigen::VectorXd rms_calc(Eigen::MatrixXd &x_kf, Eigen::MatrixXd &x_true);

double avg_rms(Eigen::VectorXd &rms_vec);