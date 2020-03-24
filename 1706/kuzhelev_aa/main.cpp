#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>
#include <iomanip>
#include <random>
#include <ctime>
#include <string>
//#include <python.h>

#define MAX_PATH 1000000000000000000

void Generate_matrix(std::vector<std::vector<int>> &matr, int size) {
	matr.resize(size);
	for (size_t i = 0; i < matr.size(); ++i) {
		matr[i].resize(size);
	}
	std::mt19937 random(static_cast<unsigned int>(time(0)));
	for (size_t i = 0; i < matr.size(); ++i) {
		for (size_t j = i + 1; j < matr.size(); ++j) {
			int num = random() % 100;
			matr[i][j] = num;
			matr[j][i] = num;
		}

	}
}

void Print_Matrix(const std::vector<std::vector<int>> &matrix) {
	
	if (matrix.size() <= 30) {
		for (const auto &row : matrix) {
			for (const auto &num : row) {
				std::cout <<std::setw(5) << num << " ";
			}
			std::cout << std::endl;
		}
	}
	else{
		std::ofstream ofs;
		ofs.open("matrix.txt");
		for (const auto &row : matrix) {
			for (const auto &num : row) {
				ofs << num << " ";
			}
			ofs << std::endl;
		}
	}
}

void Input_Matrix(std::vector<std::vector<int>> &matr, int size) {
	matr.resize(size);
	for (size_t i = 0; i < matr.size(); ++i) {
		matr[i].resize(size);
	}
	for (size_t i = 0; i < matr.size(); ++i) {
		for (size_t j = i + 1; j < matr.size(); ++j) {
			int num;
			std::cerr << "Insert distance " << i + 1 << "-" << j + 1 << ":";
			std::cin >> num;
			matr[i][j] = num;
			matr[j][i] = num;
		}

	}
}

void Print_Vector(const std::vector<int>& vec) {
	std::cout << vec[0] << ": ";
	for (size_t i = 1; i < vec[vec.size() - 1]; ++i)
		std::cout << vec[i] << " ";
	std::cout << std::endl;
}

void Read_Matrix(std::vector<std::vector<int>> &matr, size_t size) {
	matr.resize(size);
	for (size_t i = 0; i < matr.size(); ++i) {
		matr[i].resize(size, 0);
	}
	std::ifstream in;
	in.open("../../graph.txt");
	int x, y, weight;
	while (in >> x >> y >> weight) {
		matr[x][y] = matr[y][x] = weight;
	}
	in.close();
}

//void Generate_Connected_Graph(std::vector<std::vector<int>> &matr,const int& size) {
//	std::ofstream out;
//	out.open("../../graph.txt");
//	out << size;
//	out.close();
//	Py_SetProgramName(L"../../ss.py");
//	Py_Initialize();
//	PyRun_SimpleFile(fopen("../../ss.py", "r"), "ss.py");
//	matr.resize(size);
//	for (size_t i = 0; i < matr.size(); ++i) {
//		matr[i].resize(size, 0);
//	}
//	std::ifstream in;
//	in.open("../../graph.txt");
//	int x, y, weight;
//	while (in >> x >> y >> weight) {
//		matr[x][y] = matr[y][x] = weight;
//	}
//	in.close();
//}

int main() {
	int t;
	std::ofstream out("log.txt");
	std::vector<std::vector<int>> matr;
	std::cerr << "Insert peaks count:";
	std::cin >> t;
	while (t < 1) {
		std::cerr << "Peaks number can't be less than 1!\nInsert correct peaks count:";
		std::cin >> t;
	}
	int key = 1;
	if (t > 30)std::cout.rdbuf(out.rdbuf());
	std::cerr << "Input matrix by yourself [1] [default]:\nGenerate matrix(Full graph) [2] :\n\nChoose one option:";
	std::cin >> key;
	switch (key) {
	case 2:
		Generate_matrix(matr, t);
		break;
	default:
		Input_Matrix(matr, t);
	}

	std::vector<uint64_t> visited(t, 1);
	std::vector<uint64_t> min_path(t, MAX_PATH);
	Print_Matrix(matr);
	int begin_index;
	std::cerr << "Insert peak number:";
	std::cin >> begin_index;
	while (begin_index >= matr.size()) {
		std::cerr << "Index can't be more than:" << matr.size() << "\nInsert correct peak number [0, " << matr.size() - 1 << "]:";
		std::cin >> begin_index;
	}
	min_path[begin_index] = 0;
	uint64_t min_index, min;
	auto begin = std::chrono::steady_clock::now();
	do {
		min_index = MAX_PATH;
		min = MAX_PATH;
		for (size_t i = 0; i < visited.size(); ++i) {
			if (visited[i] == 1 && (min_path[i] < min)) {
				min = min_path[i];
				min_index = i;
			}
		}
		if (min_index != MAX_PATH) {
			for (size_t i = 0; i < visited.size(); ++i) {
				if (matr[min_index][i] > 0) {
					int temp = min + matr[min_index][i];
					if (temp < min_path[i]) {
						min_path[i] = temp;
					}
				}
			}
			visited[min_index] = 0;
		}
	} while (min_index < MAX_PATH);

	auto end = std::chrono::steady_clock::now();

	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

	for (size_t i = 0; i < t; ++i) {
		int end = i != begin_index ? i : i < t - 1 ? i + 1 : -1;
		if (end < 0)break;
		std::vector<int> rebuild(t + 1);
		rebuild[0] = end;
		int k = 1;
		int weight = min_path[end];
		std::vector<bool> visited(t, false);
		while (end != begin_index) {
			for (size_t j = 0; j < t; ++j) {
				int temp = weight - matr[end][j];
				if (!visited[j] && temp == min_path[j]) {
					visited[j] = true;
					weight = temp;
					end = j;
					rebuild[k] = j;
					++k;
					j = 0;
				}
			}
		}
		rebuild[t] = k;
		Print_Vector(rebuild);

	}
	std::cout << std::endl;
	for (const auto row : min_path) {
		std::cout << row << " ";
	}
	std::cout << std::endl;
	std::ofstream out1("../../info.txt");
	out1 << matr.size() << " " << begin_index;
	std::cerr << "Time spent to sequence algorithm:" << elapsed_ms.count() << " milliseconds" << std::endl;
	return 0;
}