#include <iostream>
#include <chrono>
#include <vector>
#include <fstream>
#include <iomanip>
#include <random>
#include <ctime>
#include <string>
#include <stdexcept>
//#include <python.h>

#define MAX_PATH 1000000000000000000

void Generate_matrix(std::vector<int> &matr, int size) {
	matr.resize(size*size);
	int t = size;
	std::mt19937 random(static_cast<unsigned int>(time(0)));
	for (size_t i = 0; i + 1 < matr.size(); ++i) {
		if (i%size == 0 && !i) {
			++t;
		}
		if (i%t == 0 || matr[i]) {
			continue;
		}
		matr[i] = random() % 1000;
	}
}

void Print_Matrix(const std::vector<int> &matrix) {
	
	int size = sqrt(matrix.size());
	if (matrix.size() <= 900) {
		for (int i = 0; i < matrix.size(); ++i) {
			if (i%size == 0 && i) {
				std::cout << std::endl;
			}
				std::cout << std::setw(5) << matrix[i] << " ";
		}
		std::cout << std::endl;
	}
	else{
		std::ofstream ofs;
		ofs.open("matrix.txt");
		for (const auto &row : matrix) {
				ofs << row << " ";
			
		}
		ofs << std::endl;
	}
}

void Input_Matrix(std::vector<int> &matr, int size) {
	int t = size;
	matr.resize(size*size);
	for (size_t i = 0; i < matr.size(); ++i) {
		if (i%size == 0 && !i) {
			++t;
		}
		if (i%t == 0) {
			continue;
		}
		int num ;
		std::cerr << "Insert distance " << i + 1 << "-" << i*size + 1 << ":";
		std::cin >> num;
		matr[i*size] = num;

	}
}

void Print_Vector(const std::vector<int>& vec) {
	std::cout << vec[0] << ": ";
	for (size_t i = 1; i < vec[vec.size() - 1]; ++i)
		std::cout << vec[i] << " ";
	std::cout << std::endl;
}

void Read_Matrix(std::vector<int> &matr, size_t size, const std::string &matr_info) {
	matr.resize(size*size);

	std::ifstream in;
	in.open(matr_info);
	int x, y, weight;
	while (in >> x >> y >> weight) {
		matr[x*size + y] = matr[y*size + x] = weight;
	}
	in.close();
}



//void Generate_Connected_Graph(std::vector<int> &matr,const int& size) {
//	std::ofstream out;
//	out.open("../../graph.txt");
//	out << size;
//	out.close();
//	Py_SetProgramName(L"../../ss.py");
//	Py_Initialize();
//	PyRun_SimpleFile(fopen("../../ss.py", "r"), "ss.py");
//	matr.resize(size*size);
//	std::ifstream in;
//	in.open("../../graph.txt");
//	int x, y, weight;
//	while (in >> x >> y >> weight) {
//		matr[x*size + y] = matr[y*size + x] = weight;
//	}
//	in.close();
//}

void help(const char *name, const char *message) {
	std::cout << "\n\n" << std::string(message) + 
		"\nThis is Dijkstra algorithm application.\n\n" +
		"Please provide arguments in the following format:\n\n" +
		"  $ " + name + "<matrix_size> <begin_index> <log_file_path> <save_file>\n\n" + 
		"Where <matrix_size> <begin_index> are integer numbers, " + 
		"and <log_file_path> and <save_file> are strings.\n";
}

uint64_t minDistance(std::vector<uint64_t> dist, std::vector<bool> visited) {
	uint64_t min = UINT64_MAX, min_index;

	for (size_t i = 0; i < visited.size(); ++i) {
		if (!visited[i] && dist[i] <= min) {
			min = dist[i];
			min_index = i;
		}
	}
	return min_index;
}

void PrintPath(const std::vector<int> &parent, int i) {
	if (parent[i] == -1) {
		return;
	}
	PrintPath(parent, parent[i]);
	std::cout << " " << i;
}

void PrintSolution(const std::vector<uint64_t> &min_dist, const std::vector<int> &parent, const int &begin) {
	std::cout << "Vertex"<< std::setw(10) << "Distance" << std::setw(10) << "Path\n";
	for (int i = 0; i < parent.size(); ++i) {
		if (i == begin)continue;
		std::cout << begin << " -> " << i << "\t" << min_dist[i] << "\t\t" << begin;
		PrintPath(parent, i);
		std::cout << std::endl;
	}
}

void Sequence_Alg(const std::vector<int> &matr, std::vector<uint64_t> &min_path, std::vector<int> &parent) {
	const int t = sqrt(matr.size());
	std::vector<bool> visited(t, false);

	for (int count = 0; count < t - 1; ++count) {
		int u = minDistance(min_path, visited);

		visited[u] = true;

		for (int i = 0; i < t; ++i) {
			if (!visited[i] && matr[u*t + i] && min_path[u] + matr[u*t + i] < min_path[i]) {
				parent[i] = u;
				min_path[i] = min_path[u] + matr[u*t + i];
			}
		}
	}
}

int main(int argc, char *argv[]) {
	if (argc != 5) {
		help(argv[0], "");
		return 1;
	}
	int begin_index;
	int t;
	try {
		t = std::stoi(argv[1]);
		begin_index = std::stoi(argv[2]);
	}
	catch (std::invalid_argument) {
		help(argv[0], "ERROR: Wrong matrix size value or begin index value.\n\n");
		return 1;
	}
	if (begin_index >= t) {
		help(argv[0], "ERROR: Wrong begin index value.\n\n");
		return 1;
	}
	const std::string log = std::string(argv[3]);
	const std::string info = std::string(argv[4]);
	std::ofstream out(log);
	std::vector<int> matr;
	
	
	if (t < 1) {
		help(argv[0], "ERROR: Size of matrix can't be less than 1.\\n");
		return 1;
	}
	int key = 1;
	if(t>30)std::cout.rdbuf(out.rdbuf());
	std::cerr << "Input matrix by yourself [1] [default]:\nGenerate matrix(Full graph) [2] :\n\nChoose one option:";
	std::cin >> key;
	switch (key) {
	case 2:
		Generate_matrix(matr, t);
		break;
	/*case 3:
		Generate_Connected_Graph(matr, t);
		break;*/
	default:
		Input_Matrix(matr, t);
	}

	std::vector<bool> visited(t, false);
	std::vector<uint64_t> min_path(t, MAX_PATH);
	std::vector<int> parent(t, -1);
	Print_Matrix(matr);
	
	
	
	min_path[begin_index] = 0;

	auto begin = std::chrono::steady_clock::now();
	Sequence_Alg(matr, min_path, parent);
	auto end = std::chrono::steady_clock::now();

	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

	PrintSolution(min_path, parent, begin_index);

	std::cout << std::endl;
	for (const auto row : min_path) {
		std::cout << row << " ";
	}
	std::cout << std::endl;
	std::ofstream out1(info);
	out1 << matr.size() << " " << begin_index;
	std::cerr << "Time spent to sequence algorithm:" << elapsed_ms.count() << " milliseconds" << std::endl;
	return 0;
}	Generate_matrix(matr, t);
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