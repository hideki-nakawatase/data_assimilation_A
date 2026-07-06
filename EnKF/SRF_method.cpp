#include <constant.h>
#include <fstream>
#include <Eigen/Dense>
#include <read_csv.h>
#include <lorenz96_rk4.h>
#include <random.h>
#include <jacobi.h>
#include <covariance.h>
#include <matrix_sqrt.h>
#include <rms.h>
#include <random>

using namespace std;

int main()
{
  ofstream file("EnKF_data/SRF_data.csv");

  random_device seed_gen;
  mt19937 gen(seed_gen());

  normal_distribution<double> dist_init(0.0, 1.0);

  Eigen::MatrixXd data = readCSV("observation_data/observation_data.csv");
  Eigen::MatrixXd true_data = readCSV("true_data/true_data.csv");

  Eigen::MatrixXd x_mem(N, M);
  Eigen::VectorXd ensemble_avg(N);
  Eigen::VectorXd x_f(N);
  Eigen::MatrixXd p_f(N, N);
  Eigen::MatrixXd K_gain(N, N);
  Eigen::MatrixXd K_prime_gain(N, N);
  Eigen::VectorXd x_obs(N);
  Eigen::VectorXd obs_fluctuation(N);
  Eigen::MatrixXd H(N, N);
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(N, N);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(N, N);

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

      ensemble_avg = x_mem.rowwise().mean();

      p_f.setZero();
      for (int i = 0; i < M; i++)
      {
        p_f += (x_mem.col(i) - ensemble_avg) * (x_mem.col(i) - ensemble_avg).transpose();
      }
      p_f /= M - 1;

      if (i % 10 == 0)
      {
        H = Eigen::MatrixXd::Identity(N, N);
        delete_vec = random_number_vec(step);

        for (int delete_int : delete_vec)
        {
          H(delete_int, delete_int) = 0.0;
        }

        double inf_factor = sqrt(delta_base);
        for (int j = 0; j < M; j++)
        {
          x_mem.col(j) = ensemble_avg + inf_factor * (x_mem.col(j) - ensemble_avg);
        }

        x_obs = data.row(i / 10).transpose();

        K_gain = p_f * H.transpose() * (H * p_f * H.transpose() + R).inverse();
        K_prime_gain = K_gain * (I + matrix_sqrt(R * (H * p_f * H.transpose() + R).inverse())).inverse();

        for (int j = 0; j < M; j++)
        {
          x_mem.col(j) = ensemble_avg + K_gain * (x_obs - H * ensemble_avg) + (I - K_prime_gain * H) * (x_mem.col(j) - ensemble_avg);
        }
        ensemble_avg = x_mem.rowwise().mean();
        tmp.row(i / 10) = ensemble_avg;
      }
    }

    rms_vec = rms_calc(tmp, true_data);
    result(step) = avg_rms(rms_vec);

    file << result(step) << " ";

    cout << "step " << step << endl;
  }
  return 0;
}