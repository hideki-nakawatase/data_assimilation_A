#include <random>
#include <fstream>
#include <cstdint>

using namespace std;

int N = 100000;

int main()
{
  random_device seed_gen;
  uint32_t seed = seed_gen();
  default_random_engine engine(seed);

  normal_distribution<> dist(0.0, 1.0);

  ofstream file("random.csv");

  for (int i = 0; i < N; i++)
  {
    double value = dist(engine);
    file << value << "\t\n";
  }
  return 0;
}