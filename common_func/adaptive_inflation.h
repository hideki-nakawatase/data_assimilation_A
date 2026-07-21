#pragma once
#include <Eigen/Dense>

double adaptive_inflation(Eigen::VectorXd x_obs, Eigen::VectorXd x_avg, Eigen::MatrixXd p_f, Eigen::MatrixXd R, Eigen::MatrixXd H)
{
  double inflation;
  double tmp1 = 0;
  double tmp2 = 0;
  int obs_size = x_obs.size();

  Eigen::MatrixXd tmp;

  tmp = H * p_f * H.transpose();

  Eigen::VectorXd d = x_obs - x_avg;

  for (int i = 0; i < obs_size; i++)
  {
    tmp1 += (d(i) * d(i)) * R.inverse()(i, i);
    tmp2 += tmp(i, i) * R.inverse()(i, i);
  }

  if (tmp2 < 1e-12)
    return 1.0;

  inflation = (tmp1 - obs_size) / tmp2;

  return std::max(1.0, inflation);
}