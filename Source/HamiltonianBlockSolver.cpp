/*
 * HamiltonianBlockSolver.cpp
 *
 *  Created on: Jul 30, 2016
 *      Author: Matan
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <ctime>

#include "utils.hpp"
//#include "szHamiltonian.hpp"
//#include "szBasis.hpp"
#include "genBasis.hpp"
#include "genHamiltonian.hpp"

#include "HamiltonianBlockSolver.hpp"



extern "C" void dsaupd_(int *ido, char *bmat, int *n, char *which,
	int *nev, double *tol, double *resid, int *ncv,
	double *v, int *ldv, int *iparam, int *ipntr,
	double *workd, double *workl, int *lworkl,
	int *info);

extern "C" void dseupd_(int *rvec, char *All, int *select, double *d,
	double *v, int *ldv, double *sigma,
	char *bmat, int *n, char *which, int *nev,
	double *tol, double *resid, int *ncv, double *V, //v -> V
	int *ldV, int *iparam, int *ipntr, double *workd, //ldv -> ldV
	double *workl, int *lworkl, int *ierr);


HamiltonianBlockSolver::HamiltonianBlockSolver(genBasis<double>* b, genHamiltonian<double>* h): basis(b), Hamiltonian(h) {}

void HamiltonianBlockSolver::exactSolve(double* EigenValues) {
	exactSolve(EigenValues, NULL);
}

void HamiltonianBlockSolver::exactSolve(double* EigenValues, double** EigenVectors) {
	int len = basis->getLen();

	double** HMatrix = new double*[len];
	HMatrix[0] = new double[len*len];
	for (int i = 1; i < len; i++) HMatrix[i] = HMatrix[i - 1] + len;

	Hamiltonian->toMatrix(HMatrix);
	if (EigenVectors == NULL)
		utils::diagnolize(len, HMatrix, EigenValues);
	else {
		utils::diagnolize(len, HMatrix, EigenValues,"V");
		for (int i=0; i<len; i++)
			for (int j=0; j<len; j++) EigenVectors[i][j] = HMatrix[i][j];
	}

	printEigenValues(EigenValues,5);

	delete[] HMatrix[0];
	delete[] HMatrix;
}

int HamiltonianBlockSolver::naiveLanczos(double* EigenValues, int numOfEv, int iterations, double* baseState) {
	//return naiveLanczos(EigenValues, numOfEv, iterations, baseState, 0.00000000000001);
	return naiveLanczos(EigenValues, NULL, numOfEv, iterations, baseState);
}

int HamiltonianBlockSolver::naiveLanczos(double* EigenValues, double** EigenVectors, int numOfEv, int iterations, double* baseState) {
	return naiveLanczos(EigenValues, EigenVectors, numOfEv, iterations, baseState, 0.00000000000001);
}

int HamiltonianBlockSolver::naiveLanczos(double* EigenValues, int numOfEv, int iterations, double* baseState, double precision) {
	return naiveLanczos(EigenValues, NULL, numOfEv, iterations, baseState, precision);
}
int HamiltonianBlockSolver::naiveLanczos(double* EigenValues, double** EigenVectors, int numOfEv, int iterations, double* baseState, double precision) {
	int m = iterations;
	int len = basis->getLen();
	int iterCheckConv = 10;

	// generate empty Krylov space matrix
	printf("  generating empty Krylov space matrix\n");
	double** KMatrix = new double*[m];
	KMatrix[0] = new double[m*m];
	for (int i = 1; i < m; i++) KMatrix[i] = KMatrix[i - 1] + m;

	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++) KMatrix[i][j] = 0;

	double* currEigenValues ;
	int prevN = std::max(numOfEv, iterCheckConv);
	double* prevEigenValues= new double[prevN];
	for (int i=0; i<prevN; i++) prevEigenValues[i] = 0;

	// generate initial parameters
	double* prevState = new double[len];
	double* currentState = new double[len];
	basis->copyState(baseState, currentState);
	double* tmpState = new double[len];

	double a, b2, norm, tmp;

	// populate Krylov space matrix
	printf("  populating Krylov space matrix\n");
	// get a0 and state u1
	norm = basis->getStateNorm(currentState);
	norm = norm*norm;
	Hamiltonian->apply(currentState, tmpState); // H|u_n>
	a = basis->scalarProd(currentState, tmpState) / norm;
	KMatrix[0][0] = a;

	for (int i = 0; i < len; i++) {
		tmp = currentState[i];
		currentState[i] = tmpState[i] - a*currentState[i];
		prevState[i] = tmp;
	}

	int n = 1;
	double epsilon = 0;//0.0001;
	bool converged = false;
	while (n < m && norm > epsilon && !converged) {
		b2 = 1 / norm; // 1/<u_n-1|u_n-1>
		norm = basis->getStateNorm(currentState); // sqrt(<u_n|u_n>)
		norm = norm * norm; // <u_n|u_n>
		b2 *= norm; // <u_n|u_n>/<u_n-1|u_n-1>

		Hamiltonian->apply(currentState, tmpState); // H|u_n>
		a = basis->scalarProd(currentState, tmpState) / norm; // <u_n|H|u_n>/<u_n|u_n>

		KMatrix[n][n] = a;
		KMatrix[n][n - 1] = sqrt(b2);
		KMatrix[n - 1][n] = sqrt(b2);

		for (int i = 0; i < len; i++) {
			tmp = currentState[i];
			currentState[i] = tmpState[i] - a*currentState[i] - b2*prevState[i];
			prevState[i] = tmp;
		}
		n++;

		// check convergence
		if (n % 10 == 0 && n > numOfEv) {
			currEigenValues = new double[n];
			utils::diagnolize(n, KMatrix, currEigenValues, "C");
			converged = true;
			for (int i=0; i<numOfEv; i++) {
				converged = converged & (fabs((currEigenValues[i] - prevEigenValues[i]) / (currEigenValues[i] + prevEigenValues[i]) * 2) < precision);
				/*printf("%d %d\t%.20f\t%.20f\t%.20f\n",
										n,converged,currEigenValues[i], prevEigenValues[i],
										fabs((currEigenValues[i] - prevEigenValues[i]) / (currEigenValues[i] + prevEigenValues[i]) * 2)
										);*/
			}

			delete[] prevEigenValues;
			prevN = n;
			prevEigenValues = new double[prevN];
			for (int i=0; i< prevN; i++) prevEigenValues[i] = currEigenValues[i];
			delete[] currEigenValues;
		}
		/*if (n % 10 == 0)
			printf("    %d\n", n);*/
	}
	//if (m!=n)
	//	printf("   %d  %f\n", n, norm);
	printf("  %d iterations\n",n);
	m = n;



	//printMatrix(m, KMatrix);

	// diagnolize Krylov space matrix
	printf("  diagonalizing Krylov space matrix\n");
	double* KEigenValues = new double[m];
	double** newKMatrix = new double*[m];
	newKMatrix[0] = new double[m*m];
	for (int i = 1; i < m; i++) newKMatrix[i] = newKMatrix[i - 1] + m;

	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++) newKMatrix[i][j] = KMatrix[i][j];

	int info;
	if (EigenVectors == NULL)
		info = utils::diagnolize(m, newKMatrix, KEigenValues);
	else
		info = utils::diagnolize(m, newKMatrix, KEigenValues, "V");

	if (info != 0)
		printf("  diagonalizing result: %d\n", info);
	else {
		//printf("  saving eigenvalues\n");
		//saveEigenValues(outputPath, "Lanczos", N, k, m, EigenValues);
		for (int i=0; i<numOfEv; i++) EigenValues[i] = KEigenValues[i];

		// get EigenVectors
		if (EigenVectors != NULL) {
			basis->copyState(baseState, currentState);

			norm = basis->getStateNorm(currentState);
			norm = norm*norm;
			Hamiltonian->apply(currentState, tmpState); // H|u_n>
			a = basis->scalarProd(currentState, tmpState) / norm;



			for (int i = 0; i < len; i++) {
				for(int j=0; j<numOfEv; j++) EigenVectors[j][i] = newKMatrix[j][0] * currentState[i] / sqrt(norm);

				tmp = currentState[i];
				currentState[i] = tmpState[i] - a*currentState[i];
				prevState[i] = tmp;
			}

			n = 1;
			while (n < m) {
				b2 = 1 / norm; // 1/<u_n-1|u_n-1>
				norm = basis->getStateNorm(currentState); //sqrt(<u_n|u_n>)
				norm = norm*norm; // <u_n|u_n>
				b2 *= norm; // <u_n|u_n>/<u_n-1|u_n-1>

				Hamiltonian->apply(currentState, tmpState); // H|u_n>
				a = basis->scalarProd(currentState, tmpState) / norm; // <u_n|H|u_n>/<u_n|u_n>


				for (int i = 0; i < len; i++) {
					for(int j=0; j<numOfEv; j++) EigenVectors[j][i] += newKMatrix[j][n] * currentState[i] / sqrt(norm);

					tmp = currentState[i];
					currentState[i] = tmpState[i] - a*currentState[i] - b2*prevState[i];
					prevState[i] = tmp;
				}

				n++;

				/*if (n % 10 == 0)
					printf("    %d\n", n);*/
			}
			printf("  %d iterations\n",n);
		}
		printEigenValues(EigenValues, numOfEv);

	}
	printf("\n");

	delete[] tmpState;
	delete[] prevState;
	delete[] currentState;

	delete[] KEigenValues;
	delete[] KMatrix[0];
	delete[] KMatrix;
	delete[] newKMatrix[0];
	delete[] newKMatrix;

	return m;
}

int HamiltonianBlockSolver::arpackLanczos(double* EigenValues, int numOfEv, int iterations, double* baseState) {
	return arpackLanczos(EigenValues, NULL, numOfEv, iterations, baseState);
}

int HamiltonianBlockSolver::arpackLanczos(double* EigenValues, double** EigenVectors, int numOfEv, int iterations, double* baseState) {
	return arpackLanczos(EigenValues, EigenVectors, numOfEv, iterations, baseState, 0.0, -1);
}

int HamiltonianBlockSolver::arpackLanczos(double* EigenValues, int numOfEv, int iterations, double* baseState ,double precision, int resetIters) {
	return arpackLanczos(EigenValues, NULL, numOfEv, iterations, baseState, precision, resetIters);
}

int HamiltonianBlockSolver::arpackLanczos(double* EigenValues, double** EigenVectors, int numOfEv, int iterations, double* baseState ,double precision, int resetIters) {
	printf("  running arpack Lanczos\n");
	int nev = numOfEv; // The number of values to calculate
	int n = basis->getLen();

	int ido = 0; /* Initialization of the reverse communication parameter. */

	char bmat[2] = "I";
	char which[3] = "LM";
	double tol = precision; /* Sets the tolerance; tol<=0 specifies machine precision
	 	 	 	 	 	 	 	  0.000000000000001 got max precision;
	 	 	 	 	 	 	 	  0.0000000001 got 15 significant digits
	 	 	 	 	 	 	 	 */

	double *resid = new double[n];
	basis->copyState(baseState,resid);

	int ncv; /* The largest number of basis vectors that will be used in the Implicitly Restarted
			    Arnoldi Process. Work per major iteration is proportional to N*NCV*NCV. */
	if (resetIters <= nev) ncv = 4 * nev;
	else ncv = resetIters;
	if (ncv>n) ncv = n;

	int ldv = n;
	double *v = new double[ldv*ncv];

	int *iparam;
	iparam = new int[11]; /* An array used to pass information to the
						  routines
						  about their functional modes. */
	iparam[0] = 1; // Specifies the shift strategy (1->exact)
	if (iterations < 0)
		iparam[2] = 3 * n; // Maximum number of iterations
	else
		iparam[2] = iterations;

	iparam[6] = 1; /* Sets the mode of dsaupd.
				   1 is exact shifting,
				   2 is user-supplied shifts,
				   3 is shift-invert mode,
				   4 is buckling mode,
				   5 is Cayley mode. */

	int *ipntr;
	ipntr = new int[11]; /* Indicates the locations in the work array workd
						 where the input and output vectors in the
						 callback routine are located. */

	double *workd = new double[3 * n];
	double *workl = new double[ncv*(ncv + 8)];

	int lworkl = ncv*(ncv + 8); /* Length of the workl array */

	int info = 1; /* Passes convergence information out of the iteration routine.
					 0 - random baseState
					 1 - given baseState
	 	 	 	 	 */

	int rvec = (EigenVectors != NULL); /* Specifies that eigenvectors should not be calculated */

	int *select = new int[ncv];
	double *d = new double[2 * ncv]; /* This vector will return the eigenvalues from the second routine, dseupd. */
	double sigma;
	int ierr;
	int counter = 0;
	do {
		dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
			&ncv, v, &ldv, iparam, ipntr, workd, workl,
			&lworkl, &info);

		if ((ido == 1) || (ido == -1))
			Hamiltonian->apply(workd + ipntr[0] - 1, workd + ipntr[1] - 1);

		counter++;
		/*if (counter % 10 == 0)
			printf("    %d\n", counter);*/

	} while ((ido == 1) || (ido == -1));
	printf("  %d iterations\n",counter);

	if (info<0) {
		std::cout << "Error with dsaupd, info = " << info << "\n";
		std::cout << "Check documentation in dsaupd\n\n";
	}
	else {
		dseupd_(&rvec, "All", select, d, v, &ldv, &sigma, bmat,
			&n, which, &nev, &tol, resid, &ncv, v, &ldv,
			iparam, ipntr, workd, workl, &lworkl, &ierr);

		if (ierr != 0) {
			std::cout << "Error with dseupd, info = " << ierr << "\n";
			std::cout << "Check the documentation of dseupd.\n\n";
		}
		else if (info == 1) {
			std::cout << "Maximum number of iterations reached.\n\n";
		}
		else if (info == 3) {
			std::cout << "No shifts could be applied during implicit\n";
			std::cout << "Arnoldi update, try increasing NCV.\n\n";
		}

		/* Before exiting, we copy the solution information over to
			the arrays of the calling program, then clean up the
			memory used by this routine. For some reason, when I
			don't find the eigenvectors I need to reverse the order of
			the values. */


		if (EigenVectors == NULL)
			for (int i = 0; i<nev; i++) EigenValues [i] = d[nev - 1 - i];
		else {
			for (int i = 0; i<nev; i++) {
				for (int j=0; j<n; j++) EigenVectors[i][j] = v[i * n + j];
				EigenValues [i] = d[i];
			}
		}

		delete resid;
		delete v;
		delete iparam;
		delete ipntr;
		delete workd;
		delete workl;
		delete select;
		delete d;
	}

	printEigenValues(EigenValues, nev);
	return counter;
}

void HamiltonianBlockSolver::printEigenValues(double* EigenValues, int numOfEv) {
	printf("  Eigenvalues:\n");
	for (int i=0; i< numOfEv; i++) printf("\t%.20f\n",EigenValues[i]);
}
