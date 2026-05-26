#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include "read_csv.h"
#include "lorenz96_rk4.h"
#include "covariance.h"

using namespace std;

int days = 730;
double dt = 0.005;
// double dt = 0.05;
int time_steps = 0.2 / dt * days;
int N = 40;
double F = 8;

int main()
{
    ofstream file("../KF_data/no_KF_data.csv");
    Eigen::MatrixXd p = Eigen::MatrixXd::Identity(N, N);

    Eigen::VectorXd x(N);
    vector<vector<double>> data = readCSV("../observation_data/observation_data.csv");
    for (int i = 0; i < N; i++)
    {
        x[i] = data[0][i];
        file << x[i] << " ";
    }
    file << endl;
    Eigen::VectorXd x_f(N);
    Eigen::MatrixXd p_f(N, N);
    Eigen::MatrixXd K_gain(N, N);
    Eigen::VectorXd x_obs(N);
    Eigen::MatrixXd H = Eigen::MatrixXd::Identity(N, N);
    Eigen::MatrixXd R = Eigen::MatrixXd::Identity(N, N);

    for (int i = 1; i < time_steps / 2; i++)
    {
        p = covariance_matrix(p, x);
        x = lorenz96_rk4(N, 1, dt, F, x);
        if (i % 10 == 0)
        {
            // for (int j = 0; j < N; j++)
            // {
            //     x_obs[j] = data[i / 10][j];
            // }
            // x_f = x;
            // p_f = p;
            // K_gain = p_f * H.transpose() * (H * p_f * H.transpose() + R).inverse();
            // x = x_f + K_gain * (x_obs - H * x_f);
            // p = (Eigen::MatrixXd::Identity(N, N) - K_gain * H) * p_f;
            for (int j = 0; j < N; j++)
            {
                file << x[j] << " ";
            }
            file << endl;
        }
    }
    return 0;
}