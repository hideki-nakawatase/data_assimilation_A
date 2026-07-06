#include <vector>
#include <Eigen/Dense>
#include <fstream>
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
    ofstream file("KF_data/no_KF_data.csv");
    Eigen::VectorXd x(N);
    Eigen::MatrixXd data = readCSV("observation_data/observation_data.csv");
    Eigen::MatrixXd data_true = readCSV("true_data/true_data.csv");

    for (int i = 0; i < N; i++)
    {
        x[i] = data(0, i);
    }

    for (int i = 1; i < time_steps / 2; i++)
    {
        Eigen::MatrixXd jacobi_matrix = lorenz96_jacobi_matrix(x);
        x = lorenz96_rk4(N, 1, dt, F, x);
        if (i % 10 == 0)
        {
            for (int j = 0; j < N; j++)
                file << x(j) << " ";
        }
        file << endl;
    }
    return 0;
}