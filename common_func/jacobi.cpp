#include "jacobi.h"
#include "lorenz96_rk4.h"
#include <vector>
#include <cmath>

using namespace std;

vector<vector<double>> lorenz96_jacobi_matrix(vector<double> &x)
{
  int N = 40;
  vector<vector<double>> jacobi_matrix(N, vector<double>(N, 0.0));

  double epsilon = 1e-7;

  vector<double> x_base = lorenz96_rk4(N, 1, 0.005, 8, x);

  for (int j = 0; j < N; j++)
  {
    vector<double> x_dif = x;
    x_dif[j] += epsilon;

    vector<double> x_dif_next = lorenz96_rk4(N, 1, 0.005, 8, x_dif);

    for (int i = 0; i < N; i++)
    {
      jacobi_matrix[i][j] = (x_dif_next[i] - x_base[i]) / epsilon;
    }
  }

  return jacobi_matrix;
}