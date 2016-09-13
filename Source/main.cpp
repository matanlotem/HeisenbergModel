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
#include <sstream>
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


void arpackBenchmark() {
	double Jxy = 1, Jz = 1, Hz = 0;
	bool cyclic = true;
	int N=20;

	szBasis* basis = new szBasis(N, N/2);
	szHamiltonian* Hamiltonian = new szHamiltonian(Jxy,Jz,Hz,cyclic,basis);
	HamiltonianBlockSolver<double> * solver = new HamiltonianBlockSolver<double>(basis, Hamiltonian);

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
	double Jxy = 2, Jz = 0, Hz = 0;
	bool cyclic = true;
	int N;
	int minN = 4, maxN = 20;

	szBasis* basis;
	szHamiltonian* Hamiltonian;
	HamiltonianBlockSolver<double>* solver;

	double* EigenValues;
	double* baseState;
	double* BaseEVs = new double[3];
	double EvDiff;
	int numOfEv = 2;

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
			solver = new HamiltonianBlockSolver<double>(basis, Hamiltonian);

			if (N < 12) { // exact solve
				EigenValues = new double[basis->getCombs()];
				solver->exactSolve(EigenValues);
			}
			else { // arpack solve
				EigenValues = new double[numOfEv];
				baseState = new double[basis->getCombs()];
				basis->newState(baseState, 0);
				solver->arpackLanczos(EigenValues, numOfEv, 500, baseState);
				//solver->naiveLanczos(EigenValues, numOfEv, 500, baseState);
				delete[] baseState;
			}

			// save eigenvalues
			BaseEVs[i] = EigenValues[abs(1-i)];
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

	szBasis * basis;
	szHamiltonian * Hamiltonian;
	HamiltonianBlockSolver<double> * solver;
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
		solver = new HamiltonianBlockSolver<double>(basis, Hamiltonian);
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

void exactSolve(int N) {
	int szUp;
	double Jxy=1, Jz=1, Hz=0;

	szBasis * basis;
	szHamiltonian * Hamiltonian;
	HamiltonianBlockSolver<double> * solver;
	int numOfEv = 1<<N;
	double * EigenValues = new double[numOfEv];
	int evPointer = 0;

	for (szUp = 0; szUp <= N; szUp++) {
		basis = new szBasis(N,szUp);
		Hamiltonian = new szHamiltonian(Jxy, Jz, Hz, false, basis);
		solver = new HamiltonianBlockSolver<double>(basis, Hamiltonian);
		solver->exactSolve(EigenValues + evPointer);
		evPointer += basis->getLen();
		delete solver;
		delete Hamiltonian;
		delete basis;
	}
	std::sort(EigenValues, EigenValues + numOfEv);
	for (int i=0; i<numOfEv; i++) printf("%.2f\t",EigenValues[i]);
	printf("\n");
	delete[] EigenValues;
}

void szSolve (int N) {
	std::clock_t t00,t0,t1;

	int szUp;
	double Jxy, Jz, Hz;
	bool cyclic = false;
	int numOfEv = 1;

	Jxy = 1;
	Jz = 1;
	Hz = 0;
	szUp = N/2;

	szBasis * sz;
	szHamiltonian * szH;
	HamiltonianBlockSolver<double> * realSolver;
	double* EigenValues;
	double* baseState;

	t0 = std::clock();
	t00 = t0;
	sz = new szBasis(N, szUp);
	t1 = std::clock();
	printf("Generated Basis: %.2f seconds\n", double(t1-t0) / CLOCKS_PER_SEC);
	t0 = t1;
	szH = new szHamiltonian(Jxy,Jz,Hz,cyclic,sz);
	t1 = std::clock();
	printf("Generated Hamiltonian: %.2f seconds\n", double(t1-t0) / CLOCKS_PER_SEC);
	t0 = t1;
	EigenValues = new double[numOfEv];
	baseState = new double[sz->getLen()];

	realSolver = new HamiltonianBlockSolver<double>(sz,szH);
	sz->newState(baseState,0);
	t1 = std::clock();
	printf("Generated Solver: %.2f seconds\n", double(t1-t0) / CLOCKS_PER_SEC);
	t0 = t1;
	//realSolver->naiveLanczos(EigenValues, 1, 500, baseState);
	realSolver->arpackLanczos(EigenValues, numOfEv, 500, baseState);
	t1 = std::clock();
	printf("Naive Lanczos: %.2f seconds\n", double(t1-t0) / CLOCKS_PER_SEC);
	t0 = t1;

	delete[] EigenValues;
	delete[] baseState;
	delete realSolver;
	delete szH;
	delete sz;
	printf("Total: %.2f seconds\n", double(t1-t00) / CLOCKS_PER_SEC);
}

void szTransSolve () {
	std::clock_t t00, t0,t1;

	int N, szUp, k;
	double Jxy, Jz, Hz;

	Jxy = 1;
	Jz = 1;
	Hz = 0;
	N = 22;
	szUp = N/2;

	szTransBasis * szT;
	szTransHamiltonian * szTH;
	double* EigenValues;
	std::complex<double>* baseState;
	t00 = std::clock();
	for (k=0; k<N; k++) {
		printf("\nk=%d\n==================================\n",k);
		t0 = std::clock();
		szT = new szTransBasis(N, szUp, k);
		t1 = std::clock();
		printf("Generated Basis: %.2f seconds\n", double(t1-t0) / CLOCKS_PER_SEC);
		t0 = t1;
		szTH = new szTransHamiltonian(Jxy,Jz,Hz,szT);
		t1 = std::clock();
		printf("Generated Hamiltonian: %.2f seconds\n", double(t1-t0) / CLOCKS_PER_SEC);
		t0 = t1;
		EigenValues = new double[1];
		baseState = new std::complex<double>[szT->getLen()];

		HamiltonianBlockSolver<std::complex<double> >* complexSolver;
		complexSolver = new HamiltonianBlockSolver<std::complex<double> >(szT,szTH);
		//complexSolver->exactSolve(EigenValues);
		szT->newState(baseState,0);
		t1 = std::clock();
		printf("Generated Solver: %.2f seconds\n", double(t1-t0) / CLOCKS_PER_SEC);
		t0 = t1;
		complexSolver->naiveLanczos(EigenValues, 1, 500, baseState);
		//complexSolver->arpackLanczos(EigenValues, 1, 500, baseState);
		t1 = std::clock();
		printf("Naive Lanczos: %.2f seconds\n", double(t1-t0) / CLOCKS_PER_SEC);
		t0 = t1;

		delete[] EigenValues;
		delete[] baseState;
		delete complexSolver;
		delete szTH;
		delete szT;
	}
	printf("Total: %.2f seconds\n", double(t1-t00) / CLOCKS_PER_SEC);
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
	HamiltonianBlockSolver<double> * realSolver;
	HamiltonianBlockSolver<std::complex<double> >* complexSolver;
	double * szExactEV;
	double * szNaiveLanEV;
	double * szArpackLanEV;
	double ** szEigenVectors;
	double * szBaseState;
	double * szTmpState;
	double * szTransExactEV;
	double * szTransNaiveLanEV;
	std::complex<double> * szTransBaseState;
	int szTransInd;

	int N, k;
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
				realSolver = new HamiltonianBlockSolver<double>(sz,szH);
				szBaseState = new double[sz->getLen()];
				szTmpState = new double[sz->getLen()];
				utils::logprintf(logfile, "    SUCCESS\n");

				// exact solve
				utils::logprintf(logfile, "\tRunning exact solver");
				szExactEV = new double[sz->getLen()];
				realSolver->exactSolve(szExactEV);
				utils::logprintf(logfile, "    SUCCESS\n");

				// compare first 5 eigenvalues between exact diagonalization and arpack Lanczos
				utils::logprintf(logfile, "\tRunning arpack Lanczos solver for 5 eigenvalues, no eigenvectors");
				szArpackLanEV = new double[5];
				sz->newState(szBaseState, 0);
				realSolver->arpackLanczos(szArpackLanEV,5,500,szBaseState);
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
				realSolver->arpackLanczos(szArpackLanEV,szEigenVectors,2,500,szBaseState);
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
				realSolver->naiveLanczos(szNaiveLanEV,szEigenVectors,2,500,szBaseState);
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
				szTransNaiveLanEV = new double[1];

				for (k=0; k<N; k++) {
					utils::logprintf(logfile, "\t\tGenerated szTransBasis N=%d k=%d", N, k);
					szT = new szTransBasis(N, N/2, k);
					utils::logprintf(logfile, "    SUCCESS\n");
					utils::logprintf(logfile, "\t\tGenerated szTransHamiltonian N=%d k=%d", N, k);
					szTH = new szTransHamiltonian(Jxy,Jz,Hz,szT);
					utils::logprintf(logfile, "    SUCCESS\n");

					utils::logprintf(logfile, "\tGenerating solver");
					complexSolver = new HamiltonianBlockSolver<std::complex<double> >(szT,szTH);
					utils::logprintf(logfile, "    SUCCESS\n");
					utils::logprintf(logfile, "\t\tRunning exact solver");
					complexSolver->exactSolve(szTransExactEV + szTransInd);
					utils::logprintf(logfile, "    SUCCESS\n");

					utils::logprintf(logfile, "\t\tRunning naive Lanczos solver");
					szTransBaseState = new std::complex<double>[szT->getLen()];
					szT->newState(szTransBaseState, 0);
					complexSolver->naiveLanczos(szTransNaiveLanEV,1,500,szTransBaseState);
					delete[] szTransBaseState;
					utils::logprintf(logfile, "    SUCCESS\n");
					delete complexSolver;

					// first eigenvalues
					utils::logprintf(logfile, "\t\tComparing Naive Lanczos / exact first eigenvalues");
					if (abs(szTransNaiveLanEV[0]-szTransExactEV[szTransInd]) < tol)
						utils::logprintf(logfile, "    SUCCESS\n");
					else {
						utils::logprintf(logfile, "    FAILURE\n");
						utils::logprintf(errorfile, "%s k=%d Naive Lanczos / exact first eigenvalues don't match \n",k,paramsStr);
						utils::logprintf(errorfile, "%\t%.15f\t%.15f\n",szTransNaiveLanEV[0],szTransExactEV[szTransInd]);
						errors++;
					}

					szTransInd += szT->getLen();
					delete szTH;
					delete szT;

				}
				delete[] szTransNaiveLanEV;

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

				delete realSolver;
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
	int N = 16;
	if (argc > 1) {
		std::istringstream iss(argv[1]);
		if (iss >> N) {};
	}
	//firstEnergyGap();
	//arpackBenchmark();
	//szTransTest();
	//benchmarkSolvers();
	//tester("../Logs/tester.log","../Logs/tester.err");
	//szTransSolve();
	szSolve(N);
	//exactSolve(N);

	return 0;
}
