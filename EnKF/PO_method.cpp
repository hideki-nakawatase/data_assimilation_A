#include <vector>
#include <Eigen/Dense>
#include <random>
#include <fstream>
#include <cstdint>
#include "read_csv.h"
#include "lorenz96_rk4.h"
#include "covariance.h"
#include "constant.h"
#include "jacobi.h"
#include <random.h>
#include <rms.h>

using namespace std;

int main()
{
  ofstream file("EnKF_data/PO_data.csv");

  random_device seed_gen;
  mt19937 gen(seed_gen());

  normal_distribution<double> dist_init(0.0, 1.0);
  normal_distribution<double> dist_obs(0.0, 1.0);

  Eigen::MatrixXd x_mem(N, M);
  Eigen::VectorXd ensemble_avg(N);

  Eigen::MatrixXd data = readCSV("observation_data/observation_data.csv");
  Eigen::MatrixXd data_true = readCSV("true_data/true_data.csv");

  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(N, N);

  Eigen::VectorXd x_f(N);
  Eigen::MatrixXd p(N, N);
  Eigen::MatrixXd K_gain(N, N);
  Eigen::VectorXd x_obs(N);
  Eigen::VectorXd obs_fluctuation(N);
  Eigen::MatrixXd H(N, N);

  Eigen::MatrixXd tmp(time_steps / 2 / 10, N);
  Eigen::VectorXd rms_vec(tmp.rows());
  Eigen::VectorXd result(N + 1);

  double delta_base = delta;
  for (int step = 0; step <= N; step++)
  {
    Eigen::VectorXi delete_vec(step);

    Eigen::MatrixXd H_backup = H;

    for (int i = 0; i < M; i++)
    {
      Eigen::VectorXd fluctuation(N);
      for (int j = 0; j < N; j++)
      {
        fluctuation(j) = dist_init(gen);
      }
      x_mem.col(i) = data.row(0).transpose() + fluctuation;
    }

    ensemble_avg = x_mem.rowwise().mean();
    tmp.row(0) = ensemble_avg;

    for (int i = 1; i < time_steps / 2; i++)
    {
      for (int j = 0; j < M; j++)
      {
        x_f = x_mem.col(j);
        x_mem.col(j) = lorenz96_rk4(1, x_f);
      }

      if (i % 10 == 0)
      {
        H = Eigen::MatrixXd::Identity(N, N);
        delete_vec = random_number_vec(step);

        for (int delete_int : delete_vec)
        {
          H(delete_int, delete_int) = 0.0;
        }

        x_obs = data.row(i / 10).transpose();

        ensemble_avg = x_mem.rowwise().mean();
        p.setZero();
        for (int j = 0; j < M; j++)
        {
          p += (x_mem.col(j) - ensemble_avg) * (x_mem.col(j) - ensemble_avg).transpose();
        }
        p /= M - 1;
        p *= delta_base;

        K_gain = p * H.transpose() * (H * p * H.transpose() + R).inverse();

        for (int j = 0; j < M; j++)
        {
          for (int k = 0; k < N; k++)
          {
            obs_fluctuation(k) = dist_obs(gen);
          }
          x_mem.col(j) += K_gain * (x_obs + obs_fluctuation - H * x_mem.col(j));
        }
        ensemble_avg = x_mem.rowwise().mean();
        tmp.row(i / 10) = ensemble_avg;
      }
    }
    rms_vec = rms_calc(tmp, data_true);
    result(step) = avg_rms(rms_vec);

    if (isnan(result(step)) || result(step) > 1.0)
    {
      if (delta_base > 1.5)
      {
        file << result(step) << " ";
        delta_base = delta;
        cout << "step " << step << endl;
        continue;
      }
      cout << result(step) << endl;
      step -= 1;
      delta_base += 0.1;
      H = H_backup;
      cout << "delta " << delta_base << endl;
      continue;
    }

    file << result(step) << " ";
    delta_base = delta;
    cout << "step " << step << endl;
  }
  return 0;
}