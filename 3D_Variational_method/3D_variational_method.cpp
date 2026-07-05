#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include "read_csv.h"
#include "lorenz96_rk4.h"
#include "covariance.h"
#include "constant.h"
#include "jacobi.h"

using namespace std;

int main()
{
  ofstream file("3D_Variational_method_data/data.csv");
  Eigen::MatrixXd file_output(time_steps / 10, N);
  Eigen::MatrixXd p = Eigen::MatrixXd::Identity(N, N);
  Eigen::VectorXd x(N);
  Eigen::MatrixXd x_24h(time_steps / 20, N);
  Eigen::MatrixXd x_48h(time_steps / 20, N);
  Eigen::MatrixXd data = readCSV("observation_data/observation_data.csv");
  x = data.row(0);

  Eigen::VectorXd x_f(N);
  Eigen::MatrixXd B_base(N, N);
  B_base = alpha * Eigen::MatrixXd::Identity(N, N);
  Eigen::VectorXd x_obs(N);
  Eigen::MatrixXd H = Eigen::MatrixXd::Identity(N, N);
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(N, N);

  vector<Eigen::VectorXd> analysis_history;

  // 同化ループ
  for (int i = 1; i < time_steps / 2; i++)
  {
    x = lorenz96_rk4(N, 1, dt, F, x);
    if (i % 10 == 0)
    {
      for (int j = 0; j < N; j++)
      {
        x_obs[j] = data(i / 10, j);
      }
      x_f = x;
      x = x_f + B_base * H.transpose() * (H * B_base * H.transpose() + R).inverse() * (x_obs - H * x_f);
      analysis_history.push_back(x);
    }
  }

  Eigen::MatrixXd B_raw = Eigen::MatrixXd::Zero(N, N);
  int count = 0;

  for (size_t k = 2; k < analysis_history.size(); k++)
  {
    Eigen::VectorXd x_24 = lorenz96_rk4(N, 10, dt, F, analysis_history[k - 1]);
    Eigen::VectorXd x_48 = lorenz96_rk4(N, 20, dt, F, analysis_history[k - 2]);
    Eigen::VectorXd diff = x_48 - x_24;

    B_raw += diff * diff.transpose();
    count++;
  }
  if (count > 0)
  {
    B_raw /= count;
  }

  B_raw *= 0.25;
  x = data.row(0).transpose();

  for (int i = 0; i < N; i++)
  {
    file << x[i] << " ";
  }
  file << endl;

  for (int i = 1; i < time_steps / 2; i++)
  {
    x = lorenz96_rk4(N, 1, dt, F, x);
    if (i % 10 == 0)
    {
      for (int j = 0; j < N; j++)
      {
        x_obs[j] = data(i / 10, j);
      }
      x_f = x;
      x = x_f + B_raw * H.transpose() * (H * B_raw * H.transpose() + R).inverse() * (x_obs - H * x_f);
      for (int j = 0; j < N; j++)
      {
        file << x[j] << " ";
      }
      file << endl;
    }
  }
  return 0;
}