#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include "read_csv.h"
#include "lorenz96_rk4.h"
#include "covariance.h"
#include "constant.h"
#include "random.h"
#include "rms.h"

using namespace std;

int main()
{
    ofstream file("3D_variational_method_data/random_delete_rms_avg.csv");
    ofstream file2("3D_variational_method_data/delete_index.csv");

    Eigen::VectorXd x(N);
    Eigen::MatrixXd data = readCSV("observation_data/observation_data.csv");
    Eigen::MatrixXd data_true = readCSV("true_data/true_data.csv");

    Eigen::MatrixXd tmp(time_steps / 2 / 10, N);
    Eigen::VectorXd rms_vec(tmp.rows());
    Eigen::VectorXd result(N + 1);

    Eigen::VectorXd x_f(N);
    Eigen::MatrixXd B = Eigen::MatrixXd::Identity(N, N) * alpha * alpha;
    Eigen::MatrixXd K_gain(N, N);
    Eigen::VectorXd x_obs(N);
    Eigen::MatrixXd R = Eigen::MatrixXd::Identity(N, N);

    for (int k = 0; k <= N; k++)
    {
        Eigen::VectorXi delete_vec(k);
        for (int i = 0; i < N; i++)
        {
            x[i] = data(0, i);
            tmp(0, i) = x[i];
        }

        for (int i = 1; i < time_steps / 2; i++)
        {
            x = lorenz96_rk4(N, 1, dt, F, x);
            if (i % 10 == 0)
            {
                Eigen::MatrixXd H = Eigen::MatrixXd::Identity(N, N);
                delete_vec = random_number_vec(k);
                for (int delete_int : delete_vec)
                {
                    H(delete_int, delete_int) = 0.0;
                }
                for (int j = 0; j < N; j++)
                {
                    x_obs[j] = data(i / 10, j);
                }
                x_f = x;
                K_gain = B * H.transpose() * (H * B * H.transpose() + R).inverse();
                x = x_f + K_gain * (x_obs - H * x_f);
                for (int j = 0; j < N; j++)
                    tmp(i / 10, j) = x(j);
            }
        }
        rms_vec = rms_calc(tmp, data_true);
        result(k) = avg_rms(rms_vec);

        file << result(k) << " ";
    }
    file << endl;
    return 0;
}