#include "jacobi.h"
#include "lorenz96_rk4.h"
#include <cmath>
#include <Eigen/Dense>
#include <lorenz96_rk4.h>
#include <constant.h>

using namespace std;

Eigen::MatrixXd lorenz96_jacobi_matrix(Eigen::VectorXd &x)
{
  Eigen::MatrixXd jacobi_matrix(N, N);
  Eigen::VectorXd y(N);
  jacobi_matrix.setZero();
  Eigen::VectorXd x_next(N);
  x_next = lorenz96_rk4(N, 1, dt, F, x);

  for (int i = 0; i < N; i++)
  {
    y = x;
    y(i) = x(i) + ϵ;
    y = lorenz96_rk4(N, 1, dt, F, y);
    for (int j = 0; j < N; j++)
    {
      jacobi_matrix(j, i) = (y(j) - x_next(j)) / ϵ;
    }
  }

  return jacobi_matrix;
}