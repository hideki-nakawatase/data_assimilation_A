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
#include <adaptive_inflation.h>

using namespace std;

int main()
{
    ofstream file("EnKF_data/localized_transform_data.csv");
    ofstream file2("EnKF_data/localized_transform_L.csv");

    random_device seed_gen;
    mt19937 gen(seed_gen());

    normal_distribution<double> dist_init(0.0, 1.0);

    Eigen::MatrixXd data = readCSV("observation_data/observation_data.csv");
    Eigen::MatrixXd true_data = readCSV("true_data/true_data.csv");

    Eigen::MatrixXd x_mem(N, M);
    Eigen::VectorXd ensemble_avg(N);
    Eigen::VectorXd x_f(N);
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

    Eigen::MatrixXd result_tmp(11, N + 1);
    for (int m = 0; m <= 10; m++)
    {
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
                    Eigen::MatrixXd p_f = (X_prime * X_prime.transpose()) / (M - 1);
                    Eigen::VectorXd avg_tmp = H_new * ensemble_avg;

                    double inflation = adaptive_inflation(y_new, avg_tmp, p_f, R_new, H_new);
                    inflation = std::clamp(inflation, 1.0, 1.5);
                    X_prime *= sqrt(inflation);

                    Eigen::VectorXd ensemble_avg_new(N);
                    Eigen::MatrixXd X_prime_analysis(N, M);

                    double L = 1.0 + m * 0.5;

                    for (int g = 0; g < N; g++)
                    {
                        vector<int> local_obs_indices;
                        vector<double> dist_weights;
                        for (size_t k = 0; k < keep.size(); k++)
                        {
                            int obs_gtid = keep[k];
                            double dist = abs(g - obs_gtid);
                            if (dist > N / 2)
                                dist = N - dist;

                            if (dist <= 2 * sqrt(10 / 3) * L)
                            {
                                local_obs_indices.push_back(k);
                                double rho_val = std::exp(-(dist * dist) / (2.0 * L * L));
                                dist_weights.push_back(rho_val);
                            }
                        }

                        for (int j = 0; j < M; j++)
                        {
                            x_mem.col(j) = ensemble_avg_new + X_prime_analysis.col(j);
                        }
                        int n_loc = local_obs_indices.size();

                        if (n_loc == 0)
                        {
                            ensemble_avg_new(g) = ensemble_avg(g);
                            X_prime_analysis.row(g) = X_prime.row(g);
                            continue;
                        }
                        Eigen::MatrixXd Y_prime_local(n_loc, M);
                        Eigen::MatrixXd R_local_inv = Eigen::MatrixXd::Zero(n_loc, n_loc);
                        Eigen::VectorXd innovation_local(n_loc);

                        for (int k = 0; k < n_loc; ++k)
                        {
                            int obs_idx = local_obs_indices[k];
                            Y_prime_local.row(k) = H_new.row(obs_idx) * X_prime;

                            R_local_inv(k, k) = dist_weights[k] * (1.0 / R_new(obs_idx, obs_idx));
                            innovation_local(k) = y_new(obs_idx) - H_new.row(obs_idx) * ensemble_avg;
                        }

                        Eigen::MatrixXd I_M = Eigen::MatrixXd::Identity(M, M);
                        Eigen::MatrixXd P_tilde = ((M - 1) * I_M + Y_prime_local.transpose() * R_local_inv * Y_prime_local).inverse();

                        Eigen::VectorXd w = P_tilde * Y_prime_local.transpose() * R_local_inv * innovation_local;
                        Eigen::MatrixXd T = matrix_sqrt((M - 1) * P_tilde);

                        ensemble_avg_new(g) = ensemble_avg(g) + X_prime.row(g) * w;
                        X_prime_analysis.row(g) = X_prime.row(g) * T;
                    }

                    for (int j = 0; j < M; j++)
                    {
                        x_mem.col(j) = ensemble_avg_new + X_prime_analysis.col(j);
                    }
                    ensemble_avg = ensemble_avg_new;
                    tmp.row(i / 10) = ensemble_avg;
                }
            }

            rms_vec = rms_calc(tmp, true_data);
            result_tmp(step) = avg_rms(rms_vec);
        }
    }
    for (int i = 0; i < N + 1; i++)
    {
        Eigen::Index min_row;
        file << result_tmp.col(i).minCoeff(&min_row) << " ";
        file2 << min_row << " ";
    }
    file << endl;
    file2 << endl;
    return 0;
}