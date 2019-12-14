#pragma once
#include <iostream>
#include "mpi.h"
#include <ctime>

using namespace std;

//------ ֿאנאככוכüםûי אכדמנטעל ------//
void RandomDataInitialization(double*& pMatrixGA, double*& pVectorGA, double*& pMatrixGR, double*& pVectorGR, int& Size)
{
	pMatrixGA = new double[Size * Size];
	pVectorGA = new double[Size];
	cout << "\nInput matrix:\n";
	srand(time(0));
	for (int i = 0; i < Size; i++)
	{
		for (int j = 0; j < Size; j++)
		{
			pMatrixGA[Size * i + j] = rand() % 100 + 50;
			pMatrixGA[Size * j + i] = pMatrixGA[Size * i + j];
			pMatrixGR[Size * i + j] = pMatrixGA[Size * i + j];
			pMatrixGR[Size * j + i] = pMatrixGA[Size * i + j];
		}
		pVectorGA[i] = rand() % 100 + 50;
		pVectorGR[i] = pVectorGA[i];
	}
	pPivotPos = new int[Size];
	pPivotIter = new int[Size];
	for (int i = 0; i < Size; ++i)
	{
		pPivotPos[i] = 0;
		pPivotIter[i] = -1;
	}

}
void ExactDataInitialization(double*& pVector, double*& pPrevProcResult, double*& pPrevProcD, double*& pPrevProcG,
	int* pSendInd, int Size, int RowNum, int ProcRank) {
	pPrevProcResult = new double[RowNum];
	pPrevProcD = new double[RowNum];
	pPrevProcG = new double[RowNum];
	for (int i = 0; i < RowNum; i++) {
		pPrevProcResult[i] = 0;
		pPrevProcD[i] = 0;
		pPrevProcG[i] = -pVector[pSendInd[ProcRank] / Size + i];
	}
}
void ProcessInitialization(double*& pMatrixGA, double*& pVectorGA, double*& pResultGA, double*& pMatrixGR, double*& pVectorGR, double*& pResultGR, double*& pProcRows, double*& pProcResult,
	int& Size, int& RowNum, int ProcRank, int ProcNum, double*& pPrevProcResult, double*& pPrevProcD, double*& pProcD, double*& pPrevProcG,
	double*& pProcG, double*& tmpVec, double& eps) {
	int RestRows;

	if (ProcRank == 0) {
		do {
			cout << "\n\tEnter equations count: ";
			cin >> Size;
			if (Size < ProcNum) {
				cout << "\tCount of equations must be bigger than number of processes! \n ";
			}
		} while (Size < ProcNum);
		cout << "\tEnter a valid error value:";
		cin >> eps;
	}
	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	RestRows = Size;
	for (int i = 0; i < ProcRank; i++)
		RestRows = RestRows - RestRows / (ProcNum - i);
	RowNum = RestRows / (ProcNum - ProcRank);

	pVectorGA = new double[Size];
	pResultGA = new double[Size];
	pVectorGR = new double[Size];
	pResultGR = new double[Size];
	pProcRows = new double[RowNum * Size];
	pProcResult = new double[RowNum];
	pProcD = new double[RowNum];
	pProcG = new double[RowNum];
	tmpVec = new double[RowNum];
	if (ProcRank == 0) {
		pMatrixGA = new double[Size * Size];
		pMatrixGR = new double[Size * Size];
		RandomDataInitialization(pMatrixGA, pVectorGA, pMatrixGR, pVectorGR, Size);
	}
}
void DataDistribution(double*& pMatrix, double*& pProcRows, double*& pVector, double*& pPrevProcResult, double*& pPrevProcD,
	double*& pPrevProcG, int*& pSendNum, int*& pSendInd, int Size, int RowNum, int ProcRank, int ProcNum) {
	int RestRows = Size;
	MPI_Bcast(pVector, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	pSendInd = new int[ProcNum];
	pSendNum = new int[ProcNum];
	RowNum = (Size / ProcNum);
	pSendNum[0] = RowNum * Size;
	pSendInd[0] = 0;
	for (int i = 1; i < ProcNum; i++) {
		RestRows -= RowNum;
		RowNum = RestRows / (ProcNum - i);
		pSendNum[i] = RowNum * Size;
		pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
	}

	MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcRows,
		pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

	ExactDataInitialization(pVector, pPrevProcResult, pPrevProcD, pPrevProcG, pSendInd, Size, RowNum, ProcRank);
}
void ReceiveInfoCalculation(int*& pReceiveNum, int*& pReceiveInd, int Size, int ProcNum) {
	int i;
	int RestRows = Size;
	pReceiveNum = new int[ProcNum];
	pReceiveInd = new int[ProcNum];

	pReceiveInd[0] = 0;
	pReceiveNum[0] = Size / ProcNum;
	for (i = 1; i < ProcNum; i++) {
		RestRows -= pReceiveNum[i - 1];
		pReceiveNum[i] = RestRows / (ProcNum - i);
		pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
	}
}
void ParallelMVmulCalculation(double* pProcRows, double* pVector, double* pProcResult, int Size, int RowNum) {
	int i, j;
	for (i = 0; i < RowNum; i++) {
		pProcResult[i] = 0;
		for (j = 0; j < Size; j++)
			pProcResult[i] += pProcRows[i * Size + j] * pVector[j];
	}
}
void ParallelVVsubCalculation(double* pProcLeft, double* pProcRight, double* pProcResult, int Size, int RowNum) {
	int i;
	for (i = 0; i < RowNum; i++) {
		pProcResult[i] = pProcLeft[i] - pProcRight[i];
	}
}
void ParallelVVsumCalculation(double* pProcLeft, double* pProcRight, double* pProcResult, int RowNum) {
	int i;
	for (i = 0; i < RowNum; i++) {
		pProcResult[i] = pProcLeft[i] + pProcRight[i];
	}
}
void ParallelVVmulCalculation(double* pProcLeft, double* pProcRight, double& Result, int RowNum) {
	int i;
	double ProcResult = 0;
	for (i = 0; i < RowNum; i++) {
		ProcResult += pProcLeft[i] * pProcRight[i];
	}
	MPI_Allreduce(&ProcResult, &Result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
void ParallelDVmulCalculation(double pProcLeft, double*& pProcRight, double*& ProcResult, int RowNum) {
	int i;
	for (i = 0; i < RowNum; i++) {
		ProcResult[i] = pProcLeft * pProcRight[i];
	}
}
void ParallelFirstIter(double*& pProcRows, double*& pVector, double*& pProcG, double*& pPrevProcResult, double*& tmpVec,
	int*& pReceiveNum, int*& pReceiveInd, int*& pSendInd, int ProcRank, int RowNum, int Size) {
	double* pPrevProcResultFull = new double[Size];
	MPI_Allgatherv(pPrevProcResult, pReceiveNum[ProcRank], MPI_DOUBLE, pPrevProcResultFull,
		pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);
	ParallelMVmulCalculation(pProcRows, pPrevProcResultFull, tmpVec, Size, RowNum);
	ParallelVVsubCalculation(tmpVec, pVector + pSendInd[ProcRank] / Size, pProcG, Size, RowNum);
}
void ParallelSecondIter(double*& pProcG, double*& pPrevProcG, double*& pPrevProcResult, double*& pProcD, double*& pPrevProcD, double*& tmpVec,
	int*& pReceiveNum, int*& pReceiveInd, int*& pSendInd, int ProcRank, int RowNum, int Size) {
	double parts[2];
	ParallelVVmulCalculation(pProcG, pProcG, parts[0], RowNum);
	double resparts[2];
	ParallelVVmulCalculation(pPrevProcG, pPrevProcG, parts[1], RowNum);
	ParallelDVmulCalculation((double)(parts[0] / parts[1]), pPrevProcD, tmpVec, RowNum);

	ParallelVVsubCalculation(tmpVec, pProcG, pProcD, Size, RowNum);

}
void ParallelThirdIter(double*& pProcRows, double*& pProcG, double*& pProcD, double*& tmpVec,
	int*& pReceiveNum, int*& pReceiveInd, int*& pSendInd, double& s, int ProcRank, int RowNum, int Size) {
	double above, resabove;
	double* pProcDFull = new double[Size];
	ParallelVVmulCalculation(pProcD, pProcG, above, RowNum);
	MPI_Allgatherv(pProcD, pReceiveNum[ProcRank], MPI_DOUBLE, pProcDFull,
		pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);
	ParallelMVmulCalculation(pProcRows, pProcDFull, tmpVec, Size, RowNum);
	double beyond;
	ParallelVVmulCalculation(pProcD, tmpVec, beyond, RowNum);
	s = -(double)(above / beyond);
	delete[] pProcDFull;
}
void FourthIter(double*& pProcResult, double*& pPrevProcResult, double*& pProcD, double*& tmpVec,
	int*& pReceiveNum, int*& pReceiveInd, int*& pSendInd, double& s, int ProcRank, int RowNum, int Size) {
	ParallelDVmulCalculation(s, pProcD, tmpVec, RowNum);
	ParallelVVsumCalculation(tmpVec, pPrevProcResult, pProcResult, RowNum);

}