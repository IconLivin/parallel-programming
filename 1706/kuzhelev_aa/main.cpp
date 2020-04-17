#include <omp.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <chrono>

#define MAX_PATH 1000000000000000000

void Print_Matrix(const std::vector<std::vector<int>> &matrix) {

	if (matrix.size() <= 30) {
		std::cout << std::endl;
		for (const auto &row : matrix) {
			for (const auto &num : row) {
				std::cout << std::setw(5) << num << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	else {
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

void Print_Vector(const std::vector<int>& vec) {
	std::cout << vec[0] << ": ";
	for (size_t i = 1; i < vec[vec.size() - 1]; ++i)
		std::cout << vec[i] << " ";
	std::cout << std::endl;
}

void Read_Matrix(std::vector<std::vector<int>> &matr,size_t size) {
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

int main() {
	int t;
	int begin_index;
	std::ofstream out("log.txt");
	std::ifstream in("../../info.txt");
	in >> t >> begin_index;
	if (t > 30)std::cout.rdbuf(out.rdbuf());
	std::cout << t << " " << begin_index;
	std::vector<std::vector<int>> matr;
	Read_Matrix(matr, t);
	std::vector<uint64_t> visited(t, 1);
	std::vector<uint64_t> min_path(t, MAX_PATH);
	Print_Matrix(matr);

	min_path[begin_index] = 0;
	uint64_t min_index, min, temp;
	size_t i;
	auto begin = std::chrono::steady_clock::now();
	do {
		min_index = MAX_PATH;
		min = MAX_PATH;

		#pragma omp parallel
		{
			int local_min_index = MAX_PATH;
			int local_min = MAX_PATH;

			#pragma omp for
			for (i = 0; i < visited.size(); ++i) {
				if (visited[i] == 1 && (min_path[i] < local_min)) {
					local_min = min_path[i];
					local_min_index = i;
				}
			}
			#pragma omp critical
			{
				if (local_min < min) {
					min = local_min;
					min_index = local_min_index;
				}
			}
		}
		if (min_index != MAX_PATH) {
			#pragma omp parallel for private(temp)
			for (i = 0; i < visited.size(); ++i) {
				if (matr[min_index][i] > 0) {
					temp = min + matr[min_index][i];
					if (temp < min_path[i]) {
						min_path[i] = temp;
					}
				}

				visited[min_index] = 0;
			}
		}

	} while (min_index < MAX_PATH);

	auto end = std::chrono::steady_clock::now();

	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

	std::cout << std::endl;

	for (const auto row : min_path) {
		std::cout << row << " ";
	}
	std::cerr << std::endl;
	std::cerr << "Time spent to OpenMP algorithm:" << elapsed_ms.count() << " milliseconds" << std::endl;
	return 0;
	
}