#include "random.h"
#include <random>
#include <fstream>
#include <cstdint>

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