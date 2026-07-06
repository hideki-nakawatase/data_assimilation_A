#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <cmath>

Eigen::MatrixXd matrix_sqrt(const Eigen::MatrixXd &A)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(A);

  if (solver.info() != Eigen::Success)
  {
    std::cout << "固有値失敗" << std::endl;
  }
  Eigen::VectorXd eigenvalues = solver.eigenvalues();

  Eigen::MatrixXd V = solver.eigenvectors();

  for (int i = 0; i < eigenvalues.size(); ++i)
  {
    if (eigenvalues(i) < 0.0)
      eigenvalues(i) = 0.0;
    eigenvalues(i) = std::sqrt(eigenvalues(i));
  }

  Eigen::MatrixXd sqrt_Lambda = eigenvalues.asDiagonal();

  Eigen::MatrixXd sqrt_A = V * sqrt_Lambda * V.transpose();

  return sqrt_A;
}