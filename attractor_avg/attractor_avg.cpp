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
  vector<double> X_true(N, F), X_true_cp(N);

  // attractor
  X_true[N / 2] += 0.0001;
  X_true = lorenz96_rk4(N, 1000, dt, F, X_true);
  X_true_cp = X_true;

  // goal
  X_true = lorenz96_rk4(N, time_steps, dt, F, X_true);

  // difference
  vector<double> X;
  X = X_true;
  X[N / 2] += 0.0001;

  return 0;
}