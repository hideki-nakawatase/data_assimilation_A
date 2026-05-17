#include <iostream>
#include <fstream>
#include <vector>
#include "lorenz96_rk4.h"

using namespace std;

int N = 40;
int time_steps = 100000;
double dt = 0.005;
double F = 8;

int main()
{
  ofstream ofs("../plob2_data/lyapunov_10_points.csv");

  vector<double> X_true(N, F), X_true_cp(N);
  vector<vector<double>> d(10, vector<double>(N, 0));
  vector<double> d_sum(10, 0.0), d0(10, 0.0);
  double d_sum_avg = 0.0;

  // attractor
  X_true[N / 2] += 0.0001;
  X_true = lorenz96_rk4(N, 1000, dt, F, X_true);
  X_true_cp = X_true;

  // difference
  vector<vector<double>> X(10, X_true);
  X[0][N / 2] += 0.0001;
  X[1][0] += 0.0001;
  X[1].back() += 0.0001;
  X[2][0] += 0.0001;
  X[3].back() += 0.0001;
  X[4][12] += 0.0001;
  X[5][24] += 0.0001;
  X[6][9] += 0.0001;
  X[7][19] += 0.0001;
  X[8][29] += 0.0001;
  X[9][N / 2] -= 0.0001;

  for (int k = 0; k < 10; k++)
  {
    for (int i = 0; i < N; i++)
      d[k][i] = X[k][i] - X_true[i];

    for (int i = 0; i < N; i++)
      d0[k] += d[k][i] * d[k][i];
    d0[k] = sqrt(d0[k]);
  }

  for (int j = 1; j < time_steps; j++)
  {
    vector<double> d_length(10, 0);
    for (int k = 0; k < 10; k++)
    {
      X[k] = lorenz96_rk4(N, 1, dt, F, X[k]);
    }

    X_true = lorenz96_rk4(N, 1, dt, F, X_true);
    d_sum_avg = 0.0;
    for (int k = 0; k < 10; k++)
    {
      for (int i = 0; i < N; i++)
      {
        d[k][i] = X[k][i] - X_true[i];
        d_length[k] += d[k][i] * d[k][i];
      }
      d_length[k] = sqrt(d_length[k]);
      d_sum[k] += log(d_length[k] / d0[k]);
      ofs << d_sum[k] / dt / j << " ";
      for (int i = 0; i < N; i++)
      {
        X[k][i] = X_true[i] + d0[k] * (X[k][i] - X_true[i]) / d_length[k];
      }
      d_sum_avg += d_sum[k];
    }
    ofs << d_sum_avg / 10 / dt / j << endl;
  }
  cout << d_sum_avg / 10 / dt / time_steps << endl;
  return 0;
}