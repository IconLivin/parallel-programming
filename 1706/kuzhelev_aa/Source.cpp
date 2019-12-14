#include "Gauss.h"
#include "Gradient.h"
#include "GradientParallel.h"
#include "math.h"
#include <ctime>


void ProcessTermination(double* pMatrixGA, double* pMatrixGR, double* pVectorGA, double* pVectorGR, double* pResultGA,
	double* pResultGR, double* pResultGRSeq, double* pPrevProcResult, double* pProcResult, double* pPrevResult, double* pPrevProcD,
	double* pProcD, double* pPrevD, double* pD, double* pPrevProcG, double* pProcG, double* pPrevG, double* pG, int ProcRank) {
	delete[] pPrevProcResult; delete[] pProcResult; delete[] pPrevResult; delete[] pPrevProcD;
	delete[] pProcD;  delete[] pPrevProcG; delete[] pProcG;
	if (ProcRank == 0) {
		delete[] pMatrixGR; delete[] pMatrixGA; delete[] pVectorGR; delete[] pVectorGA; delete[] pResultGA;
		delete[] pResultGR;
		delete[] pD; delete[] pG;
	}
}


void main(int argc, char* argv[]) {
	setlocale(LC_ALL, "Russian");
	double* pMatrixGA = 0; double* pMatrixGR = 0;
	double* pVectorGA = 0; double* pVectorGR = 0;
	double* pResultGA = 0; double* pResultGR = 0; double* pResultGRSeq = 0;
	double* pPrevProcResult = 0; double* pProcResult; double* pPrevResult = 0;
	double* pPrevProcD = 0; double* pProcD = 0; double* pPrevD = 0; double* pD = 0;
	double* pPrevProcG = 0; double* pProcG = 0;	double* pPrevG = 0; double* pG = 0;
	int* pSendNum; // Количество отправленных процессу элементов
	int* pSendInd; // Индекс первого среди них
	int* pReceiveNum; // Количество элементов, которые будет отправлять данный процесс
	int* pReceiveInd; // Индекс первого среди них
	double eps;//Погрешность
	double* tmpVec = 0;
	double s;//Значение, вычисляемое на шаге 3 каждой итерации
	int Size; // Размер матрицы и стобца свободных членов
	double* pProcRows;//Строки, выделенных данном процессу
	int RowNum;//Их количество
	int ProcRank, ProcNum;
	double Start, Finish, Duration, DurationPar;

	//------ Параллельная версия метода сопряженных градиентов ------//
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	ProcessInitialization(pMatrixGA, pVectorGA, pResultGA, pMatrixGR, pVectorGR, pResultGR, pProcRows, pProcResult, Size, RowNum, ProcRank, ProcNum,
		pPrevProcResult, pPrevProcD, pProcD, pPrevProcG, pProcG, tmpVec, eps);
	DataDistribution(pMatrixGR, pProcRows, pVectorGR, pPrevProcResult, pPrevProcD, pPrevProcG, pSendNum, pSendInd, Size, RowNum, ProcRank, ProcNum);
	ReceiveInfoCalculation(pReceiveNum, pReceiveInd, Size, ProcNum);
	int flag = 1;
	Start = MPI_Wtime();

	while(flag){
		double sumpogr;
		double sum = 0;
		ParallelFirstIter(pProcRows, pVectorGR, pProcG, pPrevProcResult, tmpVec, pReceiveNum, pReceiveInd, pSendInd, ProcRank, RowNum, Size);
		ParallelSecondIter(pProcG, pPrevProcG, pPrevProcResult, pProcD, pPrevProcD, tmpVec, pReceiveNum, pReceiveInd, pSendInd, ProcRank, RowNum, Size);
		ParallelThirdIter(pProcRows, pProcG, pProcD, tmpVec, pReceiveNum, pReceiveInd, pSendInd, s, ProcRank, RowNum, Size);
		FourthIter(pProcResult, pPrevProcResult, pProcD, tmpVec, pReceiveNum, pReceiveInd, pSendInd, s, ProcRank, RowNum, Size);
		for (int i = 0; i < RowNum; i++)
			sum += fabs(pProcResult[i] - pPrevProcResult[i]);
		MPI_Reduce(&sum, &sumpogr,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		if (ProcRank == 0) {
			if (sumpogr < eps) {
				flag = 0;
			}
		}
		MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		delete pPrevProcG; delete pPrevProcD; delete pPrevProcResult;
		pPrevProcG = pProcG; pPrevProcD = pProcD; pPrevProcResult = pProcResult;
		pProcG = new double[RowNum]; pProcD = new double[RowNum]; pProcResult = new double[RowNum];
	}
	MPI_Gatherv(pPrevProcResult, RowNum, MPI_DOUBLE, pResultGR,
		pReceiveNum, pReceiveInd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	Finish = MPI_Wtime();
	DurationPar = Finish - Start;
	//------------------------------------//
	if (ProcRank == 0) {
		
		Start = MPI_Wtime();
		SerialResultCalculation(pMatrixGA, pVectorGA, pResultGA, Size);
		Finish = MPI_Wtime();
		Duration = Finish - Start;
		cout << "\n\tThe result of the sequential Gauss method:\n";
		cout << "\tOperating time:" << Duration << " sec.\n";

		ProcessInitializationGradient(pVectorGR, pPrevResult, pResultGRSeq, pPrevD, pD, pPrevG, pG, Size);
		Start = MPI_Wtime();
		SoprGradSeq(pMatrixGR, pVectorGR, pPrevResult, pResultGRSeq, pPrevD, pD, pPrevG, pG, s, eps, Size);
		Finish = MPI_Wtime();
		Duration = Finish - Start;
		cout << "\n\tThe result of the sequential conjugate gradient method:\n";
		
		cout << "\tOperating time:" << Duration << " sec.\n";
		//--------------------------------------------------------------------------------------------//
		cout << "\n\tThe result of the parallel conjugate gradient method:\n";
		
		cout << "\tOperating time:" << DurationPar << " sec.\n";
		bool IsCorrect = 1;

		for (int i = 0; i < Size; i++) {
			if (abs(pResultGA[i] - pResultGR[i]) > eps)
				IsCorrect = 0;
		}
		if (IsCorrect) cout << "\n\tThe algorithm worked correctly.\n";
		else cout << "\n\tThere are errors in the algorithm.\n";
	}
	ProcessTermination(pMatrixGA, pMatrixGR, pVectorGA, pVectorGR, pResultGA, pResultGR, pResultGRSeq, pPrevProcResult, pProcResult, pPrevResult,
		pPrevProcD, pProcD, pPrevD, pD, pPrevProcG, pProcG, pPrevG, pG, ProcRank);
	MPI_Finalize();
}