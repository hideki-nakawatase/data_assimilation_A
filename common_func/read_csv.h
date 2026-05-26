#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// CSVファイルを読み込んで2次元vectorで返す関数
std::vector<std::vector<double>> readCSV(const std::string &filename);