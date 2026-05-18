#include "lorenz96_rk4.h"
#include "covariance.h"
#include <vector>

using namespace std;

vector<vector<double>> covariance_matrix(vector<double> &x)
{
  vector<vector<double>> covariance_matrix(40, vector<double>(40));
  double x_avg = 0.0, lorenz_avg = 0.0;
  for (int i = 0; i < 40; i++)
  {
    x_avg += x[i];
    lorenz_avg += lorenz96_rk4(40, 1, 0.005, 8.0, x)[i];
  }
  x_avg /= 40;
  lorenz_avg /= 40;
  return covariance_matrix;
}