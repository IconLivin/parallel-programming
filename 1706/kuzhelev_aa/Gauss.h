#pragma once
#include <iostream>

using namespace std;

//------ Метод Гаусса ------//
int* pPivotPos;
int* pPivotIter;
void Show(double* pMatrix, double* pVector, int Size)
{
	cout << "\n";
	for (int i = 0; i < Size; i++)
	{
		cout << '|';
		for (int j = 0; j < Size; j++)
		{
			cout << pMatrix[Size * i + j] << " ";
		}
		cout << '|';
		cout << '=';
		cout << '|' << pVector[i] << '|' << "\n";
	}
}
int FindPivotRow(double* pMatrix, int Size, int Iter)
{
	int PivotRow = -1;
	double MaxValue = 0;
	for (int i = 0; i < Size; i++)
	{
		if ((pPivotIter[i] == -1) && (fabs(pMatrix[i * Size + Iter]) > MaxValue))
		{
			PivotRow = i;
			MaxValue = fabs(pMatrix[i * Size + Iter]);
		}
	}
	return PivotRow;
}
void SerialColumnElimination(double* pMatrix, double* pVector, int Size, int Iter, int PivotRow)
{
	double PivotValue, PivotFactor;
	PivotValue = pMatrix[PivotRow * Size + Iter];
	for (int i = 0; i < Size; i++)
	{
		if (pPivotIter[i] == -1)
		{
			PivotFactor = pMatrix[i * Size + Iter] / PivotValue;
			for (int j = Iter; j < Size; j++)
			{
				pMatrix[i * Size + j] -= PivotFactor * pMatrix[PivotRow * Size + j];
			}
			pVector[i] -= PivotFactor * pVector[PivotRow];
		}
	}

}
void SerialGaussianElimination(double* pMatrix, double* pVector, int Size)
{
	int Iter;
	int PivotRow;
	for (Iter = 0; Iter < Size; Iter++)
	{
		PivotRow = FindPivotRow(pMatrix, Size, Iter);
		pPivotPos[Iter] = PivotRow;
		pPivotIter[PivotRow] = Iter;
		SerialColumnElimination(pMatrix, pVector, Size, Iter, PivotRow);
	}
}
void SerialBackSubstitution(double* pMatrix, double* pVector, double* pResult, int Size)
{
	int RowIndex, Row;
	for (int i = Size - 1; i >= 0; --i)
	{
		RowIndex = pPivotPos[i];
		pResult[i] = pVector[RowIndex] / pMatrix[Size * RowIndex + i];
		pMatrix[RowIndex * Size + i] = 1;
		for (int j = 0; j < i; ++j)
		{
			Row = pPivotPos[j];
			pVector[Row] -= pMatrix[Row * Size + i] * pResult[i];
			pMatrix[Row * Size + i] = 0;
		}
	}
}
void SerialResultCalculation(double* pMatrix, double* pVector, double* pResult, int Size)
{
	SerialGaussianElimination(pMatrix, pVector, Size);
	SerialBackSubstitution(pMatrix, pVector, pResult, Size);
}
