#include <tbb/tbb.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <chrono>

#define MAX_PATH UINT64_MAX

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
	std::vector<int> visited(t, 1);
	std::vector<uint64_t> min_vals(2, UINT64_MAX);
	std::vector<uint64_t> min_path(t, MAX_PATH);
	Print_Matrix(matr);

	min_path[begin_index] = 0;
	uint64_t min_index, min, temp;
	tbb::mutex mut;
	int size = tbb::task_scheduler_init::default_num_threads();
	std::cerr << "Num of threads:" << size << std::endl;
	auto begin = std::chrono::steady_clock::now();
	do {
		min_vals = tbb::parallel_reduce(
			tbb::blocked_range<uint64_t>(0, t),
			std::vector<uint64_t>(2) = {UINT64_MAX, UINT64_MAX},
			[&](const tbb::blocked_range<uint64_t>& v, std::vector<uint64_t> local_min_vals) {
				for (int i = 0; i < v.end(); ++i) {
					if (visited[i] == 1 && min_path[i] < local_min_vals[0]) {
						local_min_vals[0] = min_path[i];
						local_min_vals[1] = i;
					}
				}
				return local_min_vals;
			},
			[&](std::vector<uint64_t> x, std::vector<uint64_t> y) {
				if (x[0] < y[0]) {
					return x;
				}
				return y;
			});
		min = min_vals[0];
		min_index = min_vals[1];
		if (min_index != MAX_PATH) {
			tbb::parallel_for(
				tbb::blocked_range<uint64_t>(0,t),
				[&](const tbb::blocked_range<uint64_t>& v) {
					for (int i = 0; i < v.end(); ++i) {
						if (matr[min_index][i] > 0) {
							mut.lock();
							temp = min + matr[min_index][i];
							min_path[i] = min_path[i] > temp ? temp : min_path[i];
							mut.unlock();
						}
						visited[min_index] = 0;
					}
			});
		}

	} while (min_index < MAX_PATH);

	auto end = std::chrono::steady_clock::now();

	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

	std::cout << std::endl;

	for (const auto row : min_path) {
		std::cout << row << " ";
	}
	std::cerr << std::endl;
	std::cerr << "Time spent to TBB algorithm:" << elapsed_ms.count() << " milliseconds" << std::endl;
	return 0;

}