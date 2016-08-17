/*
 * main.cpp
 *
 *  Created on: Jul 30, 2016
 *      Author: Matan
 */
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <ctime>
#include <algorithm>

#include "utils.hpp"
#include "szHamiltonian.hpp"
#include "szTransHamiltonian.hpp"
#include "HamiltonianBlockSolver.hpp"
#include "szBasis.hpp"
#include "szTransBasis.hpp"

void exactSolve(szTransBasis * basis, szTransHamiltonian * Hamiltonian, double* EigenValues);
void exactSolve(szTransBasis * basis, szTransHamiltonian * Hamiltonian, double* EigenValues, std::complex<double>** EigenVectors);

void exactSolve(szTransBasis * basis, szTransHamiltonian * Hamiltonian, double* EigenValues) {
	exactSolve(basis, Hamiltonian, EigenValues, NULL);
}

void exactSolve(szTransBasis * basis, szTransHamiltonian * Hamiltonian, double* EigenValues, std::complex<double>** EigenVectors) {
	int len = basis->getLen();

	std::complex<double>** HMatrix = new std::complex<double>*[len];
	HMatrix[0] = new std::complex<double>[len*len];
	for (int i = 1; i < len; i++) HMatrix[i] = HMatrix[i - 1] + len;

	Hamiltonian->toMatrix(HMatrix);
	if (EigenVectors == NULL)
		utils::diagnolize(len, HMatrix, EigenValues);
	else {
		utils::diagnolize(len, HMatrix, EigenValues,"V");
		for (int i=0; i<len; i++)
			for (int j=0; j<len; j++) EigenVectors[i][j] = HMatrix[i][j];
	}

	printf("  Eigenvalues:\n");
	for (int i=0; i < std::min(5,len); i++) printf("\t%.20f\n",EigenValues[i]);

	delete[] HMatrix[0];
	delete[] HMatrix;
}

void szTransTest() {
	double Jxy, Jz, Hz;
	bool cyclic;
	int N;

	Jxy = 1;
	Jz = 1;
	Hz = 0;
	cyclic = true;
	N = 20;
	int szTot = 10;

	// regular sz basis exact diagonalization
	std::clock_t sz0 = std::clock();
	szBasis * s = new szBasis(N,szTot);
	szHamiltonian * H = new szHamiltonian(Jxy,Jz,Hz,cyclic,s);
	HamiltonianBlockSolver solver = HamiltonianBlockSolver(s,H);
	double* Eigenvalues1 = new double[s->getCombs()];
	//solver.exactSolve(Eigenvalues1);
	std::clock_t sz1 = std::clock();

	// sz translation basis exact diagonalization
	std::clock_t szT0 = std::clock();
	szTransBasis * basis;
	szTransHamiltonian * Hamiltonian;
	double* Eigenvalues2 = new double[s->getCombs()];
	int evPointer = 0;

	for (int k=0; k<N; k++) {
		basis = new szTransBasis(N,szTot,k);
		Hamiltonian = new szTransHamiltonian(Jxy,Jz,Hz,basis);
		exactSolve(basis, Hamiltonian, Eigenvalues2 + evPointer);
		evPointer += basis->getLen();
		delete Hamiltonian;
		delete basis;
	}
	std::sort(Eigenvalues2,Eigenvalues2 + s->getCombs());
	std::clock_t szT1 = std::clock();

	/*printf("\nEigenvalues:\n");
	for (int i=0; i<s->getCombs(); i++)
		printf("\t%+.20f\t%+.20f\t%d\n", Eigenvalues1[i], Eigenvalues2[i], abs(Eigenvalues1[i] -Eigenvalues2[i]) < 0.000000000000001);/**/

	printf("sz basis:\t%.2f\nszT basis:\t%.2f\n", double(sz1 - sz0) / CLOCKS_PER_SEC, double(szT1 - szT0) / CLOCKS_PER_SEC);

	delete[] Eigenvalues1;
	delete[] Eigenvalues2;
	delete H;
	delete s;
}

void arpackBenchmark() {
	double Jxy = 1, Jz = 1, Hz = 0;
	bool cyclic = true;
	int N=20;

	szBasis* basis = new szBasis(N, N/2);
	szHamiltonian* Hamiltonian = new szHamiltonian(Jxy,Jz,Hz,cyclic,basis);
	HamiltonianBlockSolver* solver = new HamiltonianBlockSolver(basis, Hamiltonian);

	double* EigenValues;
	double* baseState = new double[basis->getCombs()];
	int iterations;

	clock_t t1;
	clock_t t2;

	char buffer [50];
	sprintf(buffer,"../Output/arpackBenchmark-N_%d-Jxy_%.1f-Jz_%.1f-Hz_%.1f-cyc_%d.xls", N, Jxy, Jz, Hz, cyclic);
	FILE * file = fopen(buffer, "w");

	// headers
	fprintf(file, "numOfEv\tresetIters\tIterations\tSolveTime");
	for (int i=0; i<5; i++) fprintf(file, "\tEigenValue%d",i+1);
	fprintf(file, "\n");


	for (int numOfEv = 1; numOfEv <=5; numOfEv++) {
		EigenValues = new double[numOfEv];
		for (int resetIters = numOfEv + 2; resetIters <= std::max(numOfEv*4,30); resetIters++) {
			t1 = std::clock();
			basis->newState(baseState, 0);
			iterations = solver->arpackLanczos(EigenValues, numOfEv, 500, baseState, 0.0, resetIters);
			t2 = std::clock();


			fprintf(file, "%d\t%d\t%d\t%.2f", numOfEv, resetIters, iterations, (double(t2 - t1) / CLOCKS_PER_SEC));
			for (int i=0; i<numOfEv; i++) fprintf(file, "\t%.20f",EigenValues[i]);
			fprintf(file, "\n");

		}
		delete[] EigenValues;
	}
	fclose(file);
	delete[] baseState;
	delete solver;
	delete Hamiltonian;
	delete basis;
}

void firstEnergyGap() {
	double Jxy = 2, Jz = 2, Hz = 0;
	bool cyclic = true;
	int N;
	int numOfEv = 3;
	int minN = 3, maxN = 26;

	szBasis* basis;
	szHamiltonian* Hamiltonian;
	HamiltonianBlockSolver* solver;
	double** NEigenValues = new double*[maxN];
	NEigenValues[0] = new double[maxN*numOfEv];
	for (int i = 1; i < maxN; i++) NEigenValues[i] = NEigenValues [i - 1] + numOfEv;

	int* NIterations = new int[maxN];
	double* NHTime = new double[maxN];
	double* NSolveTime = new double[maxN];

	double* EigenValues;
	double* baseState;

	clock_t t1;
	clock_t t2;

	// Solve for different Ns
	for (N=minN; N<=maxN; N++) {
		printf("N = %d\n",N);

		// initialize system (precalc Hamiltonian)
		t1 = std::clock();
		basis = new szBasis(N, N/2);
		Hamiltonian = new szHamiltonian(Jxy,Jz,Hz,cyclic,basis);
		solver = new HamiltonianBlockSolver(basis, Hamiltonian);
		t2 = std::clock();
		NHTime[N-1] = double(t2 - t1) / CLOCKS_PER_SEC;

		t1 = t2;
		if (N < 12) { // exact solve
			EigenValues = new double[basis->getCombs()];
			solver->exactSolve(EigenValues);
			NIterations[N-1] = 0;
		}
		else { // arpack solve
			EigenValues = new double[numOfEv];
			baseState = new double[basis->getCombs()];
			basis->newState(baseState, 0);
			//NIterations[N-1] = solver->naiveLanczos(EigenValues, 2, 500, baseState);
			NIterations[N-1] = solver->arpackLanczos(EigenValues, numOfEv, 500, baseState);
			delete[] baseState;
		}

		// save eigenvalues
		for (int i=0; i<numOfEv; i++) NEigenValues[N-1][i] = EigenValues[i];
		t2 = std::clock();
		NSolveTime[N-1] = double(t2 - t1) / CLOCKS_PER_SEC;

		delete[] EigenValues;
		delete solver;
		delete Hamiltonian;
		delete basis;
	}


	// save output
	char buffer [50];
	sprintf(buffer,"../Output/minN_%d-maxN_%d-Jxy_%.1f-Jz_%.1f-Hz_%.1f-cyc_%d.xls", minN, maxN, Jxy, Jz, Hz, cyclic);
	FILE * file = fopen(buffer, "w");

	// headers
	fprintf(file, "N\tIterations\tHTime\tSolveTime\tEvDiff");
	for (int i=0; i<numOfEv; i++) fprintf(file, "\tEigenValue%d",i+1);
	fprintf(file, "\n");

	double EvDiff;
	for (N=minN; N<=maxN; N++) {
		//calc energy gap
		EvDiff = NEigenValues[N-1][0]-NEigenValues[N-1][1];
		if (fabs(EvDiff) < 0.0000001) EvDiff = NEigenValues[N-1][1]-NEigenValues[N-1][2];

		//printf("N=%d\t%d iteratrions\tEigenvalue Diff: %.5f\n",N,NIterations[N],EvDiff);
		fprintf(file, "%d\t%d\t%.2f\t%.2f\t%.5f",N, NIterations[N-1], NHTime[N-1], NSolveTime[N-1], EvDiff);
		for (int i=0; i<numOfEv; i++) fprintf(file, "\t%.20f",NEigenValues[N-1][i]);
		fprintf(file, "\n");

		//printf("\tEigenvalue:\n");
		//for (int i=0; i<numOfEv; i++) printf("\t\t%.20f\n",NEigenValues[N-1][i]);
	}
	fclose(file);

	delete[] NEigenValues[0];
	delete[] NEigenValues;
	delete[] NIterations;
	delete[] NHTime;
	delete[] NSolveTime;
}

void firstEnergyGap2() {
	double Jxy = 2, Jz = 2, Hz = 0;
	bool cyclic = true;
	int N;
	int minN = 4, maxN = 20;

	szBasis* basis;
	szHamiltonian* Hamiltonian;
	HamiltonianBlockSolver* solver;

	double* EigenValues;
	double* baseState;
	double* BaseEVs = new double[3];
	double EvDiff;

	// open file
	char buffer [50];
	sprintf(buffer,"../Output/EnergyGap/minN_%d-maxN_%d-Jxy_%.1f-Jz_%.1f-Hz_%.1f-cyc_%d.xls", minN, maxN, Jxy, Jz, Hz, cyclic);
	FILE * file = fopen(buffer, "w");

	fprintf(file, "N\tEvDiff\tSz=+1\tSz=0\tSz=-1\n");


	// Solve for different Ns
	for (N=minN; N<=maxN; N++) {
		printf("N = %d\n",N);
		for (int i=0; i<3; i++) {
			// initialize system (precalc Hamiltonian)
			basis = new szBasis(N, N/2 - i + 1);
			Hamiltonian = new szHamiltonian(Jxy,Jz,Hz,cyclic,basis);
			solver = new HamiltonianBlockSolver(basis, Hamiltonian);

			if (N < 12) { // exact solve
				EigenValues = new double[basis->getCombs()];
				solver->exactSolve(EigenValues);
			}
			else { // arpack solve
				EigenValues = new double[1];
				baseState = new double[basis->getCombs()];
				basis->newState(baseState, 0);
				solver->arpackLanczos(EigenValues, 1, 500, baseState);
				delete[] baseState;
			}

			// save eigenvalues
			BaseEVs[i] = EigenValues[0];
			delete[] EigenValues;
			delete solver;
			delete Hamiltonian;
			delete basis;
		}
		// save output
		//EvDiff = BaseEVs[1] + BaseEVs[2] - 2*BaseEVs[0];
		EvDiff = BaseEVs[2] - BaseEVs[1];
		fprintf(file, "%d\t%.5f\t%.20f\t%.20f\t%.20f\n",N,EvDiff,BaseEVs[0],BaseEVs[1],BaseEVs[2]);
	}

	fclose(file);
}

void benchmarkSolvers() {
	double Jxy, Jz, Hz;
	bool cyclic;
	int N;

	Jxy = 3;
	Jz = 2;
	Hz = 0;
	cyclic = true;

	int minN=10, maxN=22;
	int numOfEv = 1;
	int combs;

	szBasis * basis;
	szHamiltonian * Hamiltonian;
	HamiltonianBlockSolver * solver;
	double* Eigenvalues1 = new double[numOfEv];
	double* Eigenvalues2 = new double[numOfEv];
	double** EigenVectors = NULL;

	double* baseState;


	char buffer [50];
	sprintf(buffer,"../Output/benchmarkSolvers/minN_%d-maxN_%d-Jxy_%.1f-Jz_%.1f-Hz_%.1f-cyc_%d.xls", minN, maxN, Jxy, Jz, Hz, cyclic);
	FILE * file = fopen(buffer, "w");
	fprintf(file, "\tArpack Lanczos\t\t\tNaive Lanczos\t\t\n");
	fprintf(file, "N\tIterations\tTime\tEv\tIterations\tTime\tEv\n");

	std::clock_t t0,t1,t2;
	int iter1, iter2;

	for (N=minN; N<=maxN; N++) {
		printf("%d\n",N);
		basis = new szBasis(N, N/2);
		Hamiltonian = new szHamiltonian(Jxy,Jz,Hz,cyclic,basis);
		solver = new HamiltonianBlockSolver(basis, Hamiltonian);
		baseState = new double[basis->getCombs()];
		EigenVectors = new double*[numOfEv];
		EigenVectors[0] = new double[numOfEv * basis->getCombs()];
		for (int i=1; i<numOfEv; i++) EigenVectors[i] = EigenVectors[i-1] + basis->getCombs();

		t0 = std::clock();
		basis->newState(baseState,0);
		iter1 = solver->arpackLanczos(Eigenvalues1,EigenVectors,numOfEv,500,baseState);
		t1 = std::clock();
		basis->newState(baseState,0);
		iter2 = solver->naiveLanczos(Eigenvalues2,EigenVectors,numOfEv,500,baseState);
		t2 = std::clock();
		fprintf(file, "%d\t%d\t%.2f\t%f\t%d\t%.2f\t%f\n",N, iter1,double(t1 - t0) / CLOCKS_PER_SEC, Eigenvalues1[0], iter2,double(t2 - t1) / CLOCKS_PER_SEC, Eigenvalues2[0]);

		delete[] EigenVectors[0];
		delete[] EigenVectors;
		delete solver;
		delete Hamiltonian;
		delete basis;
	}
	fclose(file);
	delete[] Eigenvalues1;
	delete[] Eigenvalues2;
}

void compareSolvers() {
	double Jxy, Jz, Hz;
	bool cyclic;
	int N;

	Jxy = 1;
	Jz = 1;
	Hz = 0;
	cyclic = true;
	N = 10;

	szBasis * basis = new szBasis(N, N/2);
	szHamiltonian * Hamiltonian = new szHamiltonian(Jxy,Jz,Hz,cyclic,basis);
	//Hamiltonian->print();
	int combs = basis->getCombs();

	HamiltonianBlockSolver solver = HamiltonianBlockSolver(basis, Hamiltonian);
	int numOfEv = 4;
	double* baseState = new double[combs];
	double* tmpState = new double[combs];
	double* EigenValues;
	double** EigenVectors;

	EigenVectors= new double*[numOfEv];
	EigenVectors[0] = new double[numOfEv * combs];
	for (int i=1; i<numOfEv; i++) EigenVectors[i] = EigenVectors[i-1] + combs;

	// naive Lanczos
	basis->newState(baseState, 0);
	EigenValues = new double[numOfEv];
	solver.naiveLanczos(EigenValues, EigenVectors, numOfEv, 500, baseState);
	printf("\n");
	for (int i=0; i<numOfEv; i++) {
		for (int j=0; j<6; j++) printf("%.7f\t", EigenVectors[i][j]);
		Hamiltonian->apply(EigenVectors[i],tmpState);
		printf("\n\tnorm: %.3f\tEigenvalue: %.6f\n\n",
				basis->scalarProd(EigenVectors[i],EigenVectors[i]),
				tmpState[1]/EigenVectors[i][1]);
	}
	delete[] EigenValues;
	printf("\n");

	// arpack Lanczos
	EigenValues = new double[numOfEv];
	basis->newState(baseState, 0);
	solver.arpackLanczos(EigenValues, EigenVectors, numOfEv, 100, baseState);
	printf("\n");
	for (int i=0; i<numOfEv; i++) {
		for (int j=0; j<6; j++) printf("%.7f\t", EigenVectors[i][j]);
		Hamiltonian->apply(EigenVectors[i],tmpState);
		printf("\n\tnorm: %.3f\tEigenvalue: %.6f\n\n",
				basis->scalarProd(EigenVectors[i],EigenVectors[i]),
				tmpState[1]/EigenVectors[i][1]);
	}
	delete[] EigenValues;
	printf("\n");
	delete[] baseState;
	delete[] EigenVectors[0];
	delete[] EigenVectors;

	// exact solve
	EigenVectors = new double*[combs];
	EigenVectors[0] = new double[combs * combs];
	for (int i=1; i<combs; i++) EigenVectors[i] = EigenVectors[i-1] + combs;
	EigenValues = new double[combs];

	solver.exactSolve(EigenValues, EigenVectors);
	for (int i=0; i<numOfEv; i++) {
		for (int j=0; j<6; j++) printf("%.7f\t", EigenVectors[i][j]);
		Hamiltonian->apply(EigenVectors[i],tmpState);
		printf("\n\tnorm: %.3f\tEigenvalue: %.6f\n",
				basis->scalarProd(EigenVectors[i],EigenVectors[i]),
				tmpState[1]/EigenVectors[i][1]);
	}

	delete[] EigenValues;
	delete[] EigenVectors[0];
	delete[] EigenVectors;

	delete[] tmpState;
}

void tester(const char* logfile, const char* errorfile) {
	/*
	 * 	test odd+even N, Jz=1, Jxy=2,1, Hz=0,3
	 *
	 *	szBasis
	 *		generate basis
	 *		generate Hamiltonian
	 *	szTransBasis
	 *		generate basis
	 *		generate Hamiltonian
	 *	HamitlonianBlockSolver
	 *		exact
	 *		naiveLanczos
	 *		arpackLanczos
	 *		compare base state for all 3
	 *		compare 5 eigenvalues for arpack vs exact
	 *		validate eigenvector-eigenvalue match
	 *	TransHamiltonan
	 *		exact
	 *		compare results with szBasis exact
	 */

	double tol = 0.0000000001;
	int errors = 0;
	bool test;
	char * paramsStr = new char[10];

	utils::clearLog(logfile);
	utils::clearLog(errorfile);

	szBasis * sz;
	szTransBasis * szT;
	szHamiltonian * szH;
	szTransHamiltonian * szTH;
	HamiltonianBlockSolver * solver;
	double * szExactEV;
	double * szNaiveLanEV;
	double * szArpackLanEV;
	double ** szEigenVectors;
	double * szBaseState;
	double * szTmpState;
	double * szTransExactEV;
	int szTransInd;

	int N, szUp, k;
	double Jxy, Jz, Hz;
	bool cyclic = true;
	Jz = 1;

	for (Jxy=1; Jxy<=2; Jxy++) {
		for (Hz=0; Hz<=3; Hz+=3) {
			for (N=9; N<=10; N++) {
				sprintf(paramsStr, "(%d,%d,%d,%d)", (int) Jxy, (int) Jz, (int) Hz, N);
				utils::logprintf(logfile, "Testing Jxy=%d Jz=%d Hz=%d N=%d %s\n", (int) Jxy, (int) Jz, (int) Hz, N, paramsStr);

				utils::logprintf(logfile, "\tGenerating szBasis");
				sz = new szBasis(N, N/2);
				utils::logprintf(logfile, "    SUCCESS\n");
				utils::logprintf(logfile, "\tGenerating szHamiltonian");
				szH = new szHamiltonian(Jxy,Jz,Hz,cyclic,sz);
				utils::logprintf(logfile, "    SUCCESS\n");

				utils::logprintf(logfile, "\tGenerating solver");
				solver = new HamiltonianBlockSolver(sz,szH);
				szBaseState = new double[sz->getLen()];
				szTmpState = new double[sz->getLen()];
				utils::logprintf(logfile, "    SUCCESS\n");

				// exact solve
				utils::logprintf(logfile, "\tRunning exact solver");
				szExactEV = new double[sz->getLen()];
				solver->exactSolve(szExactEV);
				utils::logprintf(logfile, "    SUCCESS\n");

				// compare first 5 eigenvalues between exact diagonalization and arpack Lanczos
				utils::logprintf(logfile, "\tRunning arpack Lanczos solver for 5 eigenvalues, no eigenvectors");
				szArpackLanEV = new double[5];
				sz->newState(szBaseState, 0);
				solver->arpackLanczos(szArpackLanEV,5,500,szBaseState);
				utils::logprintf(logfile, "    SUCCESS\n");

				utils::logprintf(logfile, "\t\tComparing results");
				test = true;
				for (int i=0; i<5; i++) {test = test && (abs(szArpackLanEV[i]-szExactEV[i]) < tol);}
				if (test)
					utils::logprintf(logfile, "    SUCCESS\n");
				else {
					utils::logprintf(logfile, "    FAILURE\n");
					utils::logprintf(errorfile, "%s First 5 eigenvalues don't match \n",paramsStr);
					for (int i=0; i<5; i++)
						utils::logprintf(errorfile, "\t%.15f\t%.15f\n", szArpackLanEV[i], szExactEV[i]);
					errors++;
				}
				delete[] szArpackLanEV;


				// calculate arpack Lanczos and naive Lanczos first eigenvectors
				//		validate eigenvectors eigenvalues
				// 		compare first eigenvalues for arpack Lanczos, naive Lanczos and exact solve
				szEigenVectors = new double*[2];
				szEigenVectors[0] = new double[sz->getLen()*2];
				szEigenVectors[1] = szEigenVectors[0] + sz->getLen();

				// Arpack Lanczos eigenvectors eigenvalues
				utils::logprintf(logfile, "\tRunning arpack Lanczos solver for 2 eigenvalues with eigenvectors");
				szArpackLanEV = new double[2];
				sz->newState(szBaseState, 0);
				solver->arpackLanczos(szArpackLanEV,szEigenVectors,2,500,szBaseState);
				utils::logprintf(logfile, "    SUCCESS\n");
				utils::logprintf(logfile, "\t\tComparing results");
				test = true;
				for (int i=0; i<2; i++) {
					szH->apply(szEigenVectors[i], szTmpState);
					int l=0;
					while (abs(szTmpState[l]) < tol) l++;
					test = test && (abs(szTmpState[l]/szEigenVectors[i][l] - szArpackLanEV[i]) < tol);
				}
				if (test)
					utils::logprintf(logfile, "    SUCCESS\n");
				else {
					utils::logprintf(logfile, "    FAILURE\n");
					utils::logprintf(errorfile, "%s Arpack Lanczos eigenvectors-eigenvalues don't match \n",paramsStr);
					errors++;
				}

				// Naive Lanczos eigenvectors eigenvalues
				utils::logprintf(logfile, "\tRunning naive Lanczos solver for 2 eigenvalues with eigenvectors");
				szNaiveLanEV = new double[2];
				sz->newState(szBaseState, 0);
				solver->naiveLanczos(szNaiveLanEV,szEigenVectors,2,500,szBaseState);
				utils::logprintf(logfile, "    SUCCESS\n");
				utils::logprintf(logfile, "\t\tComparing results");
				test = true;
				for (int i=0; i<2; i++) {
					szH->apply(szEigenVectors[i], szTmpState);
					int l=0;
					while (abs(szTmpState[l]) < tol) l++;
					test = test && (abs(szTmpState[l]/szEigenVectors[i][l] - szNaiveLanEV[i]) < tol);
				}
				if (test)
					utils::logprintf(logfile, "    SUCCESS\n");
				else {
					utils::logprintf(logfile, "    FAILURE\n");
					utils::logprintf(errorfile, "%s Naive Lanczos eigenvectors-eigenvalues don't match \n",paramsStr);
					errors++;
				}

				// first eigenvalues
				utils::logprintf(logfile, "\tComparing Naive Lanczos / Arpack Lanczos / exact first eigenvalues");
				if ((abs(szArpackLanEV[0]-szExactEV[0]) < tol) && (abs(szNaiveLanEV[0]-szExactEV[0]) < tol))
					utils::logprintf(logfile, "    SUCCESS\n");
				else {
					utils::logprintf(logfile, "    FAILURE\n");
					utils::logprintf(errorfile, "%s Naive Lanczos / Arpack Lanczos / exact first eigenvalues don't match \n",paramsStr);
					utils::logprintf(errorfile, "%\t%.15f\t%.15f\t%.15f\n",szNaiveLanEV[0],szArpackLanEV[0],szExactEV[0]);
					errors++;
				}

				delete[] szEigenVectors[0];
				delete[] szEigenVectors;
				delete[] szArpackLanEV;
				delete[] szNaiveLanEV;


				// szTrans basis
				szTransExactEV = new double[sz->getLen()];
				szTransInd = 0;

				for (k=0; k<N; k++) {
					utils::logprintf(logfile, "\t\tGenerated szTransBasis N=%d k=%d", N, k);
					szT = new szTransBasis(N, N/2, k);
					utils::logprintf(logfile, "    SUCCESS\n");
					utils::logprintf(logfile, "\t\tGenerated szTransHamiltonian N=%d k=%d", N, k);
					szTH = new szTransHamiltonian(Jxy,Jz,Hz,szT);
					utils::logprintf(logfile, "    SUCCESS\n");

					utils::logprintf(logfile, "\t\tRunning exact solver");
					exactSolve(szT, szTH, szTransExactEV + szTransInd);
					utils::logprintf(logfile, "    SUCCESS\n");

					szTransInd += szT->getLen();
					delete szTH;
					delete szT;
				}

				utils::logprintf(logfile, "\tComparing sz basis and szTrans basis eigenvalues");
				std::sort(szTransExactEV,szTransExactEV + sz->getLen());
				test = true;
				for (int i=0; i<sz->getLen(); i++)
					test = test && (abs(szTransExactEV[i] - szExactEV[i]) < tol);
				if (test)
					utils::logprintf(logfile, "    SUCCESS\n");
				else {
					utils::logprintf(logfile, "    FAILURE\n");
					utils::logprintf(errorfile, "%s sz basis and szTrans basis eigenvalues don't match \n",paramsStr);
					errors++;
				}


				delete[] szExactEV;
				delete[] szTransExactEV;

				delete solver;
				delete szH;
				delete sz;
			}
		}
	}

	if (!errors) printf("SUCCESS!!!\n");
	else printf("FAILURE: %d errors\n", errors);

}

int main(int argc, char* argv[])
{
	//firstEnergyGap2();
	//arpackBenchmark();
	//compareSolvers();
	//szTransTest();
	//benchmarkSolvers();
	tester("../Logs/tester.log","../Logs/tester.err");

	return 0;
}
