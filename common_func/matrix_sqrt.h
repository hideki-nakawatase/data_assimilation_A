#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <cmath>

Eigen::MatrixXd matrix_sqrt(const Eigen::MatrixXd &A)
{
  if (!A.allFinite())
  {
    std::cerr << "【エラー】matrix_sqrt に NaN または Inf が渡されました！" << std::endl;
    // 処理を安全に倒すため、サイズに合わせた単位行列を返してクラッシュを防ぐ
    return Eigen::MatrixXd::Identity(A.rows(), A.cols());
  }

  Eigen::MatrixXd A_sym = 0.5 * (A + A.transpose());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(A_sym);

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