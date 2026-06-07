#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <Eigen/Dense>

Eigen::MatrixXd readCSV(const std::string &filename)
{
    Eigen::MatrixXd data;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "エラー: ファイルを開けませんでした: " << filename << std::endl;
        return data;
    }

    // 行ベクトル（std::vector）を格納する動的配列
    std::vector<std::vector<double>> matrix_data;
    std::string line;

    while (std::getline(file, line))
    {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string value;

        // スペース区切りで数値を読み込み
        while (std::getline(ss, value, ' '))
        {
            // 空の文字列をスキップ（連続するスペース対策）
            if (value.empty())
                continue;
            row.push_back(std::stod(value));
        }

        // 行が空でなければ追加
        if (!row.empty())
        {
            matrix_data.push_back(row);
        }
    }
    file.close();

    // 読み込んだデータが空の場合は空の行列を返す
    if (matrix_data.empty())
        return data;

    // 行数と列数を確定させる
    Eigen::Index rows = matrix_data.size();
    Eigen::Index cols = matrix_data[0].size();

    // Eigenの行列サイズをリサイズ
    data.resize(rows, cols);

    // std::vector から Eigen::MatrixXd へデータをコピー
    for (Eigen::Index i = 0; i < rows; ++i)
    {
        for (Eigen::Index j = 0; j < cols; ++j)
        {
            data(i, j) = matrix_data[i][j];
        }
    }

    return data;
}
