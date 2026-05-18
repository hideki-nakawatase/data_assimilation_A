#include "lorenz96_rk4.h"
#include "random.h"
#include <vector>
#include <fstream>

using namespace std;

int days = 730;
double dt = 0.005;
int time_steps = 0.2 / dt * days;
int N = 40;
double F = 8;

int main()
{
  vector<double> x(N, F), x_obs(N);
  ofstream file_true("../true_data/true_data.csv");
  ofstream file_obs("../observation_data/observation_data.csv");
  x[N / 2] += 0.0001;
  x = lorenz96_rk4(N, time_steps / 2, dt, F, x);

  for (int i = 0; i < time_steps / 2; i++)
  {
    if (i % 10 == 0)
    {
      x_obs = x;
      for (int k = 0; k < N; k++)
      {
        file_true << x[k] << " ";
        x_obs[k] += random_normal_distribution();
        file_obs << x_obs[k] << " ";
      }
      file_true << endl;
      file_obs << endl;
    }
    x = lorenz96_rk4(N, 1, dt, F, x);
  }
  return 0;
}