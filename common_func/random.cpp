#include "random.h"
#include <random>
#include <fstream>
#include <cstdint>
#include <constant.h>
#include <Eigen/Dense>

using namespace std;

double random_normal_distribution()
{
  random_device seed_gen;
  uint32_t seed = seed_gen();
  default_random_engine engine(seed);

  normal_distribution<> dist(0.0, 1.0);

  double result = dist(engine);

  return result;
}

Eigen::VectorXd random_fluctuation(double min, double max)
{
  Eigen::VectorXd result(N);
  random_device seed_gen;
  mt19937 gen(seed_gen());

  uniform_real_distribution<double> dis(min, max);

  for (int i = 0; i < N; i++)
  {
    result(i) = dis(seed_gen);
  }
  return result;
}