#include "jacobi.h"
#include "covariance.h"
#include <vector>
#include <Eigen/Dense>

using namespace std;

Eigen::MatrixXd covariance_matrix(Eigen::MatrixXd &p, Eigen::VectorXd &x)
{
  Eigen::MatrixXd jacobi_matrix = lorenz96_jacobi_matrix(x);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(40, 40);
  Eigen::MatrixXd p_next = (I + 0.005 * jacobi_matrix) * p * (I + 0.005 * jacobi_matrix).transpose();
  return p_next;
}