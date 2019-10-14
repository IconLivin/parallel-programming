#include <iostream>
#include "mpi.h"
#include <ctime>

using namespace std;

int  curr_rank_proc;
int  num_of_procs;
double  time_seq_work_alg = 0;
double  time_pp_work_alg = 0;


double* Create_vector_and_init(int size_row, int size_column)
{
	if (size_row < 1 || size_column < 1)
		return NULL;

	double* vec;
	vec = new double[size_row * size_column];
	srand(time(NULL));
	if (size_row > 4 || size_column > 4)
	{
		for (int i = 0; i < size_row * size_column; i++)vec[i] = (double)(rand() % 100) - 50.0 + 0.5;;
	}
	else
	{
		for (int i = 0; i < size_row * size_column; i++)
		{
			cout << "Input" << i << " double element: " << endl;
			cin >> vec[i];
		}
	}
	return vec;
}

void Show_vec(double* vec, int size_matr_row, int size_matr_column)
{
	if (vec != NULL)
		for (int i = 0; i < size_matr_row * size_matr_column; i++)
		{
			cout << vec[i] << " ";
			cout << endl;
		}
}
int main(int argc, char* argv[])
{
	int size_row = 5, size_column = 5;
	double* matrix_as_vector = NULL;

	double sum_el_seq = 0;
	double sum_el_pp = 0;

	double end_time_of_seq_alg = 0;
	double start_time_of_seq_alg = 0;
	double end_time_of_pp_alg = 0;
	double start_time_of_pp_alg = 0;


	double partial_summ = 0, temp_sum = 0;
	int size_el = 1;
	int size_work_of_proc = 0;

	MPI_Status stat;


	/* Начало MPI кода */

	MPI_Init(&argc, &argv);


	MPI_Comm_size(MPI_COMM_WORLD, &num_of_procs);

	MPI_Comm_rank(MPI_COMM_WORLD, &curr_rank_proc);

	if (num_of_procs < 1)
	{
		cout << "Incorrect number of processes (at least 1 must be)" << endl;
		return 0;
	}

	if (curr_rank_proc == 0)
	{
		cout << "Input size of row: " << endl;
		cin >> size_row;
		cout << "Input size of column: " << endl;
		cin >> size_column;

		size_el = size_row * size_column;
		matrix_as_vector = Create_vector_and_init(size_row, size_column);

		if (matrix_as_vector == NULL)
		{
			cout << "Incorrect input data, try again" << endl;
			return 0;
		}

		if (size_row * size_column < 1000)
		{
			cout << "Current matrix:" << endl;
			Show_vec(matrix_as_vector, size_row, size_column);
		}
		cout << endl;

		/* Подсчет суммы всех элементов матрицы (последовательная версия алгоритма): */

		start_time_of_seq_alg = MPI_Wtime();

		for (int i = 0; i < size_row * size_column; i++)
			sum_el_seq += matrix_as_vector[i];

		end_time_of_seq_alg = MPI_Wtime();
		time_seq_work_alg = end_time_of_seq_alg - start_time_of_seq_alg;

		cout << "Sum of all elements in matrix is  " << sum_el_seq << " " << endl;
		cout << "Spent time on the implementation of this algorithm (Sequence version)" << time_seq_work_alg << " ms " << endl;
		cout << endl;
		cout << "Num of procs: " << num_of_procs << endl;



		start_time_of_pp_alg = MPI_Wtime();

		size_work_of_proc = size_el / num_of_procs;

		MPI_Bcast(&size_el, 1, MPI_INT, 0, MPI_COMM_WORLD);
		for (int i = 1; i < num_of_procs; i++)
			MPI_Send(matrix_as_vector + size_work_of_proc * (i - 1), size_work_of_proc, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);


		cout << "Process with rang " << curr_rank_proc << " start own job" << endl;

		for (int i = size_work_of_proc * (num_of_procs - 1); i < size_el; i++)//подсчитывает последнюю оставшуюся часть, если однопроцесорная система, считает все сама
			temp_sum += matrix_as_vector[i];

		sum_el_pp = temp_sum;

		MPI_Reduce(&partial_summ, &temp_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		sum_el_pp += temp_sum;
		

		end_time_of_pp_alg = MPI_Wtime();

		time_pp_work_alg = end_time_of_pp_alg - start_time_of_pp_alg;
		cout << "Sum of all elements in matrix is (Parallel version): " << sum_el_pp << endl;
		cout << "Spent time on the implementation of this algorithm (Parallel version)" << time_pp_work_alg << " ms" << endl;

		if (sum_el_pp == sum_el_seq)
			cout << "Results of parallel and sequence versions are identical! " << endl;
		else
			cout << "Results of parallel and sequence versions are not identical! " << endl;

		if (time_pp_work_alg <= time_seq_work_alg)
			cout << "Parallel version faster, then sequence" << endl;
		else
			cout << "Sequence version faster, then parallel" << endl;

		delete[] matrix_as_vector;
	}
	else
	{
		double* recv_matrix_as_vector;

		MPI_Bcast(&size_el, 1, MPI_INT, 0, MPI_COMM_WORLD);

		size_work_of_proc = size_el / num_of_procs;
		recv_matrix_as_vector = new double[size_work_of_proc];
		MPI_Recv(recv_matrix_as_vector, size_work_of_proc, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &stat);

		cout << "Process with rang " << curr_rank_proc << " start own job" << endl;

		for (int i = 0; i < size_work_of_proc; i++)
			partial_summ += recv_matrix_as_vector[i];

		MPI_Reduce(&partial_summ, &temp_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		delete[] recv_matrix_as_vector;
	}
	MPI_Finalize();

	/* Конец MPI кода */
	return 0;
}