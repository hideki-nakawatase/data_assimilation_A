#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

ofstream ofs("lyapunov.csv");

int time_steps = 100000;
int N = 40;
double dt = 0.005;
double F = 8;

vector<double> lorenz96(const vector<double> &X)
{
  vector<double> dX(N);
  for (int i = 0; i < N; i++)
  {
    int p1 = (i + 1) % N;
    int m1 = (i - 1 + N) % N;
    int m2 = (i - 2 + N) % N;

    dX[i] = (X[p1] - X[m2]) * X[m1] - X[i] + F;
  }
  return dX;
}

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
  vector<double> d(N, 0), X(N, F), X_true(N, F);
  vector<double> k1(N), k2(N), k3(N), k4(N);
  vector<double> dX1(N), dX2(N), dX3(N), dX4(N), tmp(N);
  vector<double> d_ls(0);

  X_true[N / 2] += 0.01;

  for (int j = 0; j < 10000; j++)
  {
    dX1 = lorenz96(X_true);
    for (int i = 0; i < N; i++)
      k1[i] = dt * dX1[i];

    for (int i = 0; i < N; i++)
      tmp[i] = X_true[i] + k1[i] / 2;

    dX2 = lorenz96(tmp);
    for (int i = 0; i < N; i++)
      k2[i] = dt * dX2[i];

    for (int i = 0; i < N; i++)
      tmp[i] = X_true[i] + k2[i] / 2;

    dX3 = lorenz96(tmp);
    for (int i = 0; i < N; i++)
      k3[i] = dt * dX3[i];

    for (int i = 0; i < N; i++)
      tmp[i] = X_true[i] + k3[i];

    dX4 = lorenz96(tmp);
    for (int i = 0; i < N; i++)
      k4[i] = dt * dX4[i];

    for (int i = 0; i < N; i++)
      X_true[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
  }

  X = X_true;
  X[N / 2] += 0.001;

  for (int i = 0; i < N; i++)
    d[i] = X[i] - X_true[i];

  double d0 = 0;
  for (int i = 0; i < N; i++)
    d0 += d[i] * d[i];
  d0 = sqrt(d0);

  for (int j = 0; j < time_steps; j++)
  {
    // X
    dX1 = lorenz96(X);
    for (int i = 0; i < N; i++)
      k1[i] = dt * dX1[i];

    for (int i = 0; i < N; i++)
      tmp[i] = X[i] + k1[i] / 2;

    dX2 = lorenz96(tmp);
    for (int i = 0; i < N; i++)
      k2[i] = dt * dX2[i];

    for (int i = 0; i < N; i++)
      tmp[i] = X[i] + k2[i] / 2;

    dX3 = lorenz96(tmp);
    for (int i = 0; i < N; i++)
      k3[i] = dt * dX3[i];

    for (int i = 0; i < N; i++)
      tmp[i] = X[i] + k3[i];

    dX4 = lorenz96(tmp);
    for (int i = 0; i < N; i++)
      k4[i] = dt * dX4[i];

    for (int i = 0; i < N; i++)
      X[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
    // X_fix
    dX1 = lorenz96(X_true);
    for (int i = 0; i < N; i++)
      k1[i] = dt * dX1[i];

    for (int i = 0; i < N; i++)
      tmp[i] = X_true[i] + k1[i] / 2;

    dX2 = lorenz96(tmp);
    for (int i = 0; i < N; i++)
      k2[i] = dt * dX2[i];

    for (int i = 0; i < N; i++)
      tmp[i] = X_true[i] + k2[i] / 2;

    dX3 = lorenz96(tmp);
    for (int i = 0; i < N; i++)
      k3[i] = dt * dX3[i];

    for (int i = 0; i < N; i++)
      tmp[i] = X_true[i] + k3[i];

    dX4 = lorenz96(tmp);
    for (int i = 0; i < N; i++)
      k4[i] = dt * dX4[i];

    for (int i = 0; i < N; i++)
      X_true[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;

    // d
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