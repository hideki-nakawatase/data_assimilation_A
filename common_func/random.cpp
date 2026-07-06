#include "random.h"
#include <random>
#include <fstream>
#include <Eigen/Dense>
#include <cstdint>
#include <constant.h>

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

int uniform_random_number(int min, int max)
{
  random_device seed_gen;
  mt19937 engine(seed_gen());
  uniform_int_distribution<int> dist(min, max);
  return dist(engine);
}

Eigen::VectorXi random_number_vec(int index)
{
  Eigen::VectorXi delete_vec(index);
  Eigen::VectorXi rest_index = Eigen::VectorXi::LinSpaced(N, 0, N - 1);
  for (int k = 0; k < index; k++)
  {
    int delete_pos = uniform_random_number(0, N - k - 1);
    int delete_int = rest_index(delete_pos);
    delete_vec(k) = delete_int;
    rest_index(delete_pos) = rest_index(N - k - 1);
  }

  return delete_vec;
}