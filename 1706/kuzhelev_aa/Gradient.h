#pragma once
//------ Последовательный алгоритм ------//
double* MVmul(double* pMatrix, double* pVector, const int& Size) {
	double* res = new double[Size];
	for (int i = 0; i < Size; ++i) {
		res[i] = 0;
		for (int j = 0; j < Size; j++)
			res[i] += pMatrix[Size * i + j] * pVector[j];
	}
	return res;
}
double* VVsub(double* pVectorLeft, double* pVectorRight, const int& Size) {
	double* res = new double[Size];
	for (int i = 0; i < Size; ++i) {
		res[i] = pVectorLeft[i] - pVectorRight[i];
	}
	return res;
}
double* VVsum(double* pVectorLeft, double* pVectorRight, const int& Size) {
	double* res = new double[Size];
	for (int i = 0; i < Size; ++i) {
		res[i] = pVectorLeft[i] + pVectorRight[i];
	}
	return res;
}
double ScalMul(double* pVectorLeft, double* pVectorRight, const int& Size) {
	double res = 0;
	for (int i = 0; i < Size; ++i) {
		res += pVectorLeft[i] * pVectorRight[i];
	}
	return res;
}
double* VDmul(double* pVector, double pDouble, const int& Size) {
	double* res = new double[Size];
	for (int i = 0; i < Size; ++i) {
		res[i] = pVector[i] * pDouble;
	}
	return res;
}
void ProcessInitializationGradient(double*& pVector, double*& pPrevResult, double*& pResult, double*& pPrevD, double*& pD,
	double*& pPrevG, double*& pG, int Size) {
	pPrevResult = new double[Size];
	pResult = new double[Size];
	pPrevD = new double[Size];
	pD = new double[Size];
	pPrevG = new double[Size];
	pG = new double[Size];
	for (int i = 0; i < Size; i++) {
		pPrevResult[i] = 0;
		pPrevD[i] = 0;
		pPrevG[i] = -pVector[i];
	}
}
void SoprGradSeq(double*& pMatrixGrad, double*& pVectorGrad, double*& pPrevResult, double*& pResult, double*& pPrevD, double*& pD,
	double*& pPrevG, double*& pG, double& s, double eps, int Size) {
	ProcessInitializationGradient(pVectorGrad, pPrevResult, pResult, pPrevD, pD, pPrevG, pG, Size);
	int flag = 1;
	while (flag) {
		//first
		double* pLeft = MVmul(pMatrixGrad, pPrevResult, Size);
		pG = VVsub(pLeft, pVectorGrad, Size);
		delete pLeft;
		//second
		double above = ScalMul(pG, pG, Size);
		double beyond = ScalMul(pPrevG, pPrevG, Size);
		double* pRight = VDmul(pPrevD, (double)(above / beyond), Size);
		pD = VVsub(pRight, pG, Size);
		delete pRight;
		//third
		above = ScalMul(pD, pG, Size);
		double* beyondr = MVmul(pMatrixGrad, pD, Size);
		beyond = ScalMul(pD, beyondr, Size);
		delete beyondr;
		s = -(double)(above / beyond);
		//fourth
		pRight = VDmul(pD, s, Size);
		pResult = VVsum(pPrevResult, pRight, Size);
		double sum = 0;
		for (int i = 0; i < Size; i++)
			sum += fabs(pResult[i] - pPrevResult[i]);
		delete pPrevG; delete pPrevD; delete pPrevResult;
		pPrevG = pG; pPrevD = pD; pPrevResult = pResult;
		if (sum < eps) flag = 0;
	}
}
