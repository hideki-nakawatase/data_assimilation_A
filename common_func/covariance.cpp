#include "jacobi.h"
#include "covariance.h"
#include "constant.h"
#include <vector>
#include <Eigen/Dense>

using namespace std;

Eigen::MatrixXd covariance_matrix(Eigen::MatrixXd &p, Eigen::VectorXd &x, Eigen::MatrixXd &M)
{
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(40, 40);
  Eigen::MatrixXd p_next = delta * (I + 0.005 * M) * p * (I + 0.005 * M).transpose();
  return p_next;
}