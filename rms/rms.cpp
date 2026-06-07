#include "rms.h"
#include "read_csv.h"
#include "constant.h"
#include <Eigen/Dense>
#include <fstream>

using namespace std;

int main()
{
    ofstream file("KF_data/KF_rms_delta.csv");
    Eigen::MatrixXd kf_data = readCSV("KF_data/KF_data_delta.csv");
    Eigen::MatrixXd true_data = readCSV("true_data/true_data.csv");

    Eigen::VectorXd rms = rms_calc(kf_data, true_data);

    for (int i = 0; i < rms.rows(); i++)
    {
        file << i * dt * 10;
        for (int j = 0; j < rms.cols(); j++)
        {
            file << "," << rms(i, j);
        }
        file << endl;
    }
    return 0;
}