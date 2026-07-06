#include "lorenz96_rk4.h"
#include "lorenz96.h"
#include <Eigen/Dense>
#include <constant.h>

using namespace std;

Eigen::VectorXd lorenz96_rk4(int steps, Eigen::VectorXd &X)
{
  Eigen::VectorXd k1(N), k2(N), k3(N), k4(N);
  Eigen::VectorXd dX1(N), dX2(N), dX3(N), dX4(N), tmp(N), X_next(N);
  X_next = X;
  for (int j = 0; j < steps; j++)
  {
    dX1 = lorenz96(X_next, N, F);
    for (int i = 0; i < N; i++)
      k1[i] = dt * dX1[i];

    for (int i = 0; i < N; i++)
      tmp[i] = X[i] + k1[i] / 2;

    dX2 = lorenz96(tmp, N, F);
    for (int i = 0; i < N; i++)
      k2[i] = dt * dX2[i];

    for (int i = 0; i < N; i++)
      tmp[i] = X[i] + k2[i] / 2;

    dX3 = lorenz96(tmp, N, F);
    for (int i = 0; i < N; i++)
      k3[i] = dt * dX3[i];

    for (int i = 0; i < N; i++)
      tmp[i] = X[i] + k3[i];

    dX4 = lorenz96(tmp, N, F);
    for (int i = 0; i < N; i++)
      k4[i] = dt * dX4[i];

    for (int i = 0; i < N; i++)
      X_next[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
  }
  return X_next;
}