/*
 * HamiltonianBlockSolver.hpp
 *
 *  Created on: Jul 30, 2016
 *      Author: Matan
 */

#ifndef HAMILTONIANBLOCKSOLVER_HPP_
#define HAMILTONIANBLOCKSOLVER_HPP_
#include <complex>
#include "genBasis.hpp"
#include "genHamiltonian.hpp"

template <class T>
class HamiltonianBlockSolver {
private:
	genBasis<T> * basis;
	genHamiltonian<T> * Hamiltonian;
public:
	HamiltonianBlockSolver(genBasis<T> * b, genHamiltonian<T> * h);

	void exactSolve(double* EigenValues);
	void exactSolve(double* EigenValues, T** EigenVectors);

	int naiveLanczos(double* EigenValues, int numOfEv, int iterations, T* baseState);
	int naiveLanczos(double* EigenValues, T** EigenVectors, int numOfEv, int iterations, T* baseState);
	int naiveLanczos(double* EigenValues, int numOfEv, int iterations, T* baseState, double precision);
	int naiveLanczos(double* EigenValues, T** EigenVectors, int numOfEv, int iterations, T* baseState, double precision);

	int arpackLanczos(double* EigenValues, int numOfEv, int iterations, T* baseState);
	int arpackLanczos(double* EigenValues, T** EigenVectors, int numOfEv, int iterations, T* baseState);
	int arpackLanczos(double* EigenValues, int numOfEv, int iterations, T* baseState, double precision, int resetIters);
	int arpackLanczos(double* EigenValues, T** EigenVectors, int numOfEv, int iterations, T* baseState, double precision, int resetIters);

	void printEigenValues(double* EigenValues, int numOfEv);
};


template class HamiltonianBlockSolver<double>;
template class HamiltonianBlockSolver<std::complex<double> >;

#endif /* HAMILTONIANBLOCKSOLVER_HPP_ */
