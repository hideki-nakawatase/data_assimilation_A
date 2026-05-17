#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include "lorenz96_rk4.h"

using namespace std;

int time_steps = 100000;
int N = 40;
double dt = 0.005;
double F = 8;

double sum(const vector<double> &X)
{
  double sum = 0;
  for (int i = 0; i < X.size(); i++)
  {
    sum += X[i];
  }
  return sum;
}

int main()
{
  ofstream ofs("plob2_data/lyapunov.csv");

  vector<double> d(N, 0), X(N, F), X_true(N, F);
  vector<double> k1(N), k2(N), k3(N), k4(N);
  vector<double> dX1(N), dX2(N), dX3(N), dX4(N), tmp(N);
  vector<double> d_ls(0);

  X_true[N / 2] += 0.01;

  X_true = lorenz96_rk4(N, 10000, dt, F, X_true);

  X = X_true;
  X[N / 2] += 0.001;

  for (int i = 0; i < N; i++)
    d[i] = X[i] - X_true[i];

  double d0 = 0;
  for (int i = 0; i < N; i++)
    d0 += d[i] * d[i];
  d0 = sqrt(d0);

  X = lorenz96_rk4(N, time_steps, dt, F, X);

  X_true = lorenz96_rk4(N, time_steps, dt, F, X_true);

  // d
  for (int j = 0; j < time_steps; j++)
  {
    double d_length = 0;
    for (int i = 0; i < N; i++)
    {
      d[i] = X[i] - X_true[i];
      d_length += d[i] * d[i];
    }
    d_length = sqrt(d_length);
    d_ls.resize(d_ls.size() + 1);
    d_ls[j] = log(d_length / d0);
    ofs << sum(d_ls) / dt / j << endl;

    for (int i = 0; i < N; i++)
      X[i] = X_true[i] + d0 * (X[i] - X_true[i]) / d_length;
  }
  cout << "λ = " << sum(d_ls) / dt / time_steps << endl;
  return 0;
}