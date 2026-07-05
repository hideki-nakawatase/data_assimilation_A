#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include "read_csv.h"
#include "lorenz96_rk4.h"
#include "covariance.h"
#include "constant.h"
#include "jacobi.h"
#include <random.h>

using namespace std;

int main()
{
  ofstream file("3D_Variational_method_data/ensemble_data.csv");
  Eigen::MatrixXd file_output(time_steps / 10, N);
  Eigen::VectorXd x(N);
  Eigen::MatrixXd data = readCSV("observation_data/observation_data.csv");
  x = data.row(0);

  Eigen::VectorXd x_f(N);
  Eigen::MatrixXd B_base(N, N);
  Eigen::MatrixXd B_raw(N, N);
  B_base = alpha * Eigen::MatrixXd::Identity(N, N);
  Eigen::VectorXd x_obs(N);
  Eigen::MatrixXd H = Eigen::MatrixXd::Identity(N, N);
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(N, N);

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