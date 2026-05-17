#include "lorenz96.h"

using namespace std;

vector<double> lorenz96(vector<double> &X, int N, double F)
{
    {
        vector<double> dX(N);
        for (int i = 0; i < N; i++)
        {
            int p1 = (i + 1) % N;
            int m1 = (i - 1 + N) % N;
            int m2 = (i - 2 + N) % N;

            dX[i] = (X[p1] - X[m2]) * X[m1] - X[i] + F;
        }
        return dX;
    }
}
