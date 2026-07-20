#include <constant.h>
#include <fstream>
#include <Eigen/Dense>
#include <read_csv.h>
#include <lorenz96_rk4.h>
#include <random.h>
#include <jacobi.h>
#include <covariance.h>
#include <matrix_sqrt.h>
#include <rms.h>
#include <random>

using namespace std;

int main()
{
    ofstream file("EnKF_data/transform_data.csv");

    random_device seed_gen;
    mt19937 gen(seed_gen());

    normal_distribution<double> dist_init(0.0, 1.0);

    Eigen::MatrixXd data = readCSV("observation_data/observation_data.csv");
    Eigen::MatrixXd true_data = readCSV("true_data/true_data.csv");

    Eigen::MatrixXd x_mem(N, M);
    Eigen::VectorXd ensemble_avg(N);
    Eigen::VectorXd x_f(N);
    Eigen::MatrixXd p_f(N, N);
    Eigen::MatrixXd K_gain(N, N);
    Eigen::MatrixXd K_prime_gain(N, N);
    Eigen::VectorXd x_obs(N);
    Eigen::VectorXd obs_fluctuation(N);
    Eigen::MatrixXd H(N, N);
    Eigen::MatrixXd R = Eigen::MatrixXd::Identity(N, N);
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(N, N);

    Eigen::MatrixXd tmp(time_steps / 2 / 10, N);
    Eigen::VectorXd rms_vec(tmp.rows());
    Eigen::VectorXd result(N + 1);

    double delta_base = delta;
    for (int step = 0; step <= N; step++)
    {
        Eigen::VectorXi delete_vec(step);

        Eigen::MatrixXd H_backup = H;

        for (int i = 0; i < M; i++)
        {
            Eigen::VectorXd fluctuation(N);
            for (int j = 0; j < N; j++)
            {
                fluctuation(j) = dist_init(gen);
            }
            x_mem.col(i) = data.row(0).transpose() + fluctuation;
        }

        ensemble_avg = x_mem.rowwise().mean();
        tmp.row(0) = ensemble_avg;

        for (int i = 1; i < time_steps / 2; i++)
        {
            for (int j = 0; j < M; j++)
            {
                x_f = x_mem.col(j);
                x_mem.col(j) = lorenz96_rk4(1, x_f);
            }

            ensemble_avg = x_mem.rowwise().mean();

            Eigen::MatrixXd X_prime(N, M);
            for (int j = 0; j < M; j++)
            {
                X_prime.col(j) = x_mem.col(j) - ensemble_avg;
            }
            X_prime *= sqrt(delta_base);

            p_f.setZero();
            for (int i = 0; i < M; i++)
            {
                p_f += X_prime.col(i) * X_prime.col(i).transpose();
            }
            p_f /= M - 1;

            if (i % 10 == 0)
            {
                H = Eigen::MatrixXd::Identity(N, N);
                delete_vec = random_number_vec(step);
                std::vector<int> keep;

                for (int i = 0; i < H.rows(); ++i)
                {
                    auto it = std::find(delete_vec.begin(), delete_vec.end(), i);

                    if (it == delete_vec.end())
                    {
                        keep.push_back(i);
                    }
                }

                Eigen::MatrixXd H_new(keep.size(), H.cols());

                for (size_t i = 0; i < keep.size(); ++i)
                {
                    H_new.row(i) = H.row(keep[i]);
                }

                H = H_new;
                Eigen::MatrixXd R_new = Eigen::MatrixXd::Identity(keep.size(), keep.size());

                Eigen::VectorXd y_full = data.row(i / 10).transpose();
                Eigen::VectorXd y_new(keep.size());
                for (size_t k = 0; k < keep.size(); ++k)
                {
                    y_new(k) = y_full[keep[k]];
                }

                Eigen::MatrixXd X_prime(N, M);
                for (int j = 0; j < M; j++)
                {
                    X_prime.col(j) = x_mem.col(j) - ensemble_avg;
                }

                Eigen::MatrixXd Y_prime = H_new * X_prime;

                Eigen::MatrixXd I_M = Eigen::MatrixXd::Identity(M, M);
                Eigen::MatrixXd I_obs = Eigen::MatrixXd::Identity(keep.size(), keep.size());

                Eigen::MatrixXd P_tilde = (I_M - Y_prime.transpose() * (Y_prime * Y_prime.transpose() + (M - 1) * I_obs).inverse() * Y_prime) / (M - 1);

                Eigen::VectorXd innovation = y_new - H_new * ensemble_avg;
                Eigen::VectorXd w = P_tilde * Y_prime.transpose() * innovation;
                Eigen::VectorXd ensemble_avg_new = ensemble_avg + X_prime * w;

                Eigen::MatrixXd T = matrix_sqrt((M - 1) * P_tilde);
                Eigen::MatrixXd X_prime_analysis = X_prime * T;

                for (int j = 0; j < M; j++)
                {
                    x_mem.col(j) = ensemble_avg_new + X_prime_analysis.col(j);
                }

                ensemble_avg = ensemble_avg_new;
                tmp.row(i / 10) = ensemble_avg;
            }
        }

        rms_vec = rms_calc(tmp, true_data);
        result(step) = avg_rms(rms_vec);

        file << result(step) << " ";

        cout << "step " << step << endl;
    }
    return 0;
}