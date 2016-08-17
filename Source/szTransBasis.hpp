/*
 * szTransBasis.hpp
 *
 *  Created on: Jul 30, 2016
 *      Author: Matan
 */

#ifndef SZTRANSBASIS_HPP_
#define SZTRANSBASIS_HPP_

#include <complex>

#include "szBasis.hpp"

class szTransBasis {
private:
	int N;
	int SzUp;
	int K;
	int combs;
	int len;
	szBasis * szB;
	int* MapSzBasisIndex2TBasisIndex;
	std::complex<double>* MapSzBasisIndex2TNormFactor;
	int* MapTBasisIndex2SzBasisState;
	double* MapTBasisIndex2TNormFactor;
public:
	szTransBasis(int n, int szUp, int k);
	~szTransBasis();
	int getCombs();
	int getLen();
	int getN();
	int getK();
	int getSzUp();

	int szBasisIndex2TBasisIndex(int szBasisIndex);
	std::complex<double> szBasisIndex2TNormFactor(int szBasisIndex);
	int TBasisIndex2szBasisState(int TBasisIndex);
	int szBasisState2szBasisIndex(int szBasisState);
	int szBasisIndex2szBasisState(int szBasisIndex);
	double TBasisIndex2TNormFactor(int TBasisIndex);

	void newState(std::complex<double>* state);
	void newState(std::complex<double>* state, int TBasisIndex);
	void newState(std::complex<double>* state, int TBasisIndex, std::complex<double> phase);
	void copyState(std::complex<double>* inputState, std::complex<double>* outputState);
	std::complex<double> scalarProd(std::complex<double>* state1, std::complex<double>* state2);
	double getStateNorm(std::complex<double>* state);
	double normalizeState(std::complex<double>* state);
	void printState(std::complex<double>* state);
};



#endif /* SZTRANSBASIS_HPP_ */
