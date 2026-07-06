#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <cmath>
#include "read_csv.h"
#include "lorenz96_rk4.h"
#include "covariance.h"
#include "constant.h"
#include "jacobi.h"
#include "rms.h"
#include "random.h"

using namespace std;

int main()
{
    ofstream file("KF_data/KF_data_delete.csv");
    if (!file)
        cout << "file not found" << endl;

    Eigen::MatrixXd p = Eigen::MatrixXd::Identity(N, N);
    Eigen::VectorXd x(N);

    Eigen::MatrixXd data = readCSV("observation_data/observation_data.csv");
    Eigen::MatrixXd data_true = readCSV("true_data/true_data.csv");

    Eigen::MatrixXd tmp(time_steps / 2 / 10, N);
    Eigen::VectorXd rms_vec(tmp.rows());
    Eigen::VectorXd result(N + 1);

    Eigen::VectorXd x_f(N);
    Eigen::MatrixXd p_f(N, N);
    Eigen::MatrixXd K_gain(N, N);
    Eigen::VectorXd x_obs(N);
    Eigen::MatrixXd H = Eigen::MatrixXd::Identity(N, N);
    Eigen::MatrixXd R = Eigen::MatrixXd::Identity(N, N);

    double delta_base = 1.01;
    for (int step = 0; step <= N; step++)
    {
        Eigen::VectorXi delete_vec(step);

        Eigen::MatrixXd H_backup = H;

        for (int i = 0; i < N; i++)
        {
            x[i] = data(0, i);
            tmp(0, i) = x[i];
        }
        p = Eigen::MatrixXd::Identity(N, N);

        for (int i = 1; i < time_steps / 2; i++)
        {
            Eigen::MatrixXd jacobi_matrix = lorenz96_jacobi_matrix(x);
            p = covariance_matrix(p, jacobi_matrix);
            p = 0.5 * (p + p.transpose());
            x = lorenz96_rk4(N, 1, dt, F, x);

            if (i % 10 == 0)
            {
                H = Eigen::MatrixXd::Identity(N, N);
                delete_vec = random_number_vec(step);
                for (int delete_int : delete_vec)
                {
                    H(delete_int, delete_int) = 0.0;
                }
                for (int j = 0; j < N; j++)
                {
                    x_obs[j] = data(i / 10, j);
                }
                x_f = x;
                p_f = p * delta_base;
                K_gain = p_f * H.transpose() * (H * p_f * H.transpose() + R).inverse();
                x = x_f + K_gain * (x_obs - H * x_f);
                p = (Eigen::MatrixXd::Identity(N, N) - K_gain * H) * p_f;
                for (int j = 0; j < N; j++)
                    tmp(i / 10, j) = x(j);
            }
        }
        rms_vec = rms_calc(tmp, data_true);
        result(step) = avg_rms(rms_vec);

        if (isnan(result(step)) || result(step) > 1.5)
        {
            if (delta_base > 2.0)
            {
                file << result(step) << " ";
                delta_base = 1.01;
                cout << "step " << step << endl;
                continue;
            }
            cout << result(step) << endl;
            step -= 1;
            delta_base += 0.1;
            H = H_backup;
            cout << "delta " << delta_base << endl;
            continue;
        }

        file << result(step) << " ";
        delta_base = 1.01;
        cout << "step " << step << endl;
    }
    return 0;
}