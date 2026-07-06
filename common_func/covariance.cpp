#include "jacobi.h"
#include "covariance.h"
#include "constant.h"
#include <vector>
#include <Eigen/Dense>

using namespace std;

Eigen::MatrixXd covariance_matrix(Eigen::MatrixXd &p, Eigen::MatrixXd &M)
{
  Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(N, N);
  Eigen::MatrixXd p_next = M * p * M.transpose();
  return p_next;
}