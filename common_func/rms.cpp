#include "rms.h"
#include "constant.h"
#include "read_csv.h"
#include <Eigen/Dense>

using namespace std;

Eigen::VectorXd rms_calc(Eigen::MatrixXd &x_kf, Eigen::MatrixXd &x_true)
{
    int kf_row = x_kf.rows();
    int kf_col = x_kf.cols();
    Eigen::VectorXd rms_vec(kf_row);

    for (int i = 0; i < kf_row; i++)
    {
        double sum = 0;
        for (int j = 0; j < kf_col; j++)
        {
            sum += (x_kf(i, j) - x_true(i, j)) * (x_kf(i, j) - x_true(i, j));
        }
        rms_vec(i) = pow(sum / kf_col, 0.5);
    }

    return rms_vec;
}

double avg_rms(Eigen::VectorXd &rms_vec)
{
    double avg = 0.0;
    int spin_up_lim = 60 / 0.2 * 4 / 10;
    for (int i = spin_up_lim; i < rms_vec.rows(); i++)
        avg += rms_vec[i];
    avg /= rms_vec.rows() - spin_up_lim;

    return avg;
}