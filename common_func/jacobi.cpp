#include "jacobi.h"
#include "lorenz96_rk4.h"
#include <cmath>
#include <Eigen/Dense>

using namespace std;

Eigen::MatrixXd lorenz96_jacobi_matrix(Eigen::VectorXd &x)
{
  int N = 40;
  Eigen::MatrixXd jacobi_matrix(N, N);
  jacobi_matrix.setZero();

  for (int i = 0; i < N; i++)
  {
    int k = (i - 1 + N) % N;
    int l = (i - 2 + N) % N;
    int m = (i + 1) % N;

    jacobi_matrix(i, i) = -1.0;
    jacobi_matrix(i, k) = x[m] - x[l];
    jacobi_matrix(i, l) = -x[k];
    jacobi_matrix(i, m) = x[k];
  }

  return jacobi_matrix;
}