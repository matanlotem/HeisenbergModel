/*
 * HamiltonianBlockSolver.hpp
 *
 *  Created on: Jul 30, 2016
 *      Author: Matan
 */

#ifndef HAMILTONIANBLOCKSOLVER_HPP_
#define HAMILTONIANBLOCKSOLVER_HPP_

class HamiltonianBlockSolver {
private:
	szBasis * basis;
	szHamiltonian * Hamiltonian;
public:
	HamiltonianBlockSolver(szBasis * b, szHamiltonian * h);

	void exactSolve(double* EigenValues);
	void exactSolve(double* EigenValues, double** EigenVectors);

	int naiveLanczos(double* EigenValues, int numOfEv, int iterations, double* baseState);
	int naiveLanczos(double* EigenValues, double** EigenVectors, int numOfEv, int iterations, double* baseState);
	int naiveLanczos(double* EigenValues, int numOfEv, int iterations, double* baseState, double precision);
	int naiveLanczos(double* EigenValues, double** EigenVectors, int numOfEv, int iterations, double* baseState, double precision);

	int arpackLanczos(double* EigenValues, int numOfEv, int iterations, double* baseState);
	int arpackLanczos(double* EigenValues, double** EigenVectors, int numOfEv, int iterations, double* baseState);
	int arpackLanczos(double* EigenValues, int numOfEv, int iterations, double* baseState, double precision, int resetIters);
	int arpackLanczos(double* EigenValues, double** EigenVectors, int numOfEv, int iterations, double* baseState, double precision, int resetIters);

	void printEigenValues(double* EigenValues, int numOfEv);
};


#endif /* HAMILTONIANBLOCKSOLVER_HPP_ */
