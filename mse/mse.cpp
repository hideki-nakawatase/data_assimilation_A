#include <iostream>
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
    ofstream ofs("../plob2_data/diff_distance.csv");

    vector<double> X_true(N, F), d_mse_avg_ls(N);
    vector<vector<double>> d(10, vector<double>(N, 0));
    vector<double> d0(10, 0.0), d_mse(10, 0.0), d_mse_cp(10, 0.0);

    // attractor
    X_true[N / 2] += 0.0001;
    X_true = lorenz96_rk4(N, 1000, dt, F, X_true);

    // difference
    vector<vector<double>> X(10, X_true);
    X[0][N / 2] += 0.0001;
    X[1][0] -= 0.0001;
    X[2][0] += 0.0001;
    X[3].back() += 0.0001;
    X[4][12] += 0.0001;
    X[5][24] += 0.0001;
    X[6][9] += 0.0001;
    X[7][19] += 0.0001;
    X[8][29] += 0.0001;
    X[9][N / 2] -= 0.0001;

    for (int k = 0; k < 10; k++)
    {
        for (int i = 0; i < N; i++)
            d[k][i] = X[k][i] - X_true[i];

        for (int i = 0; i < N; i++)
            d0[k] += d[k][i] * d[k][i];
        d0[k] = sqrt(d0[k]);
    }

    for (int j = 1; j < time_steps; j++)
    {
        vector<double> d_length(10, 0);
        double d_mse_avg = 0.0;
        double d_mse_avg_cp;
        for (int k = 0; k < 10; k++)
        {
            X[k] = lorenz96_rk4(N, 1, dt, F, X[k]);
        }

        X_true = lorenz96_rk4(N, 1, dt, F, X_true);
        for (int k = 0; k < 10; k++)
        {
            for (int i = 0; i < N; i++)
            {
                d[k][i] = X[k][i] - X_true[i];
                d_length[k] += d[k][i] * d[k][i];
            }
            d_length[k] = sqrt(d_length[k]);
            d_mse[k] = log(d_length[k] / sqrt(N) / d0[k]);
            ofs << d_mse[k] << " ";
            d_mse_avg += d_mse[k];
        }
        ofs << d_mse_avg / 10 << endl;
        d_mse_avg_ls[j] = d_mse_avg - d_mse_avg_cp;
        d_mse_avg_cp = d_mse_avg;
    }
    return 0;
}