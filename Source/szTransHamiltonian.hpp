/*
 * szTransHamiltonian.hpp
 *
 *  Created on: Jul 30, 2016
 *      Author: Matan
 */

#ifndef SZTRANSHAMILTONIAN_HPP_
#define SZTRANSHAMILTONIAN_HPP_
#include <string>
#include "szTransBasis.hpp"

class szTransHamiltonian {
protected:
	int N;
	int K;
	int SzUp;
	int len;
	double Jxy;
	double Jz;
	double Hz;
	szTransBasis * basis;
private:
	int** ladderMatrix;
	std::complex<double>** ladderValue;
	std::complex<double>* stateValues;

	int szLadderUp(int i, int basisState);
	int szLadderDown(int i, int basisState);
	int szLadderUpDown(int i, int basisState);
	int szLadderDownUp(int i, int basisState);
	double szSpinValue(int i, int basisState);
public:
	szTransHamiltonian(double jxy, double jz, double hz, szTransBasis * b);
	~szTransHamiltonian();

	void ladderUpDown(int basisState, int* outputStates, std::complex<double>* outputValues);
	void ladderDownUp(int basisState, int* outputStates, std::complex<double>* outputValues);
	std::complex<double> spinValue(int basisState);
	std::complex<double> pairSpinValue(int basisState);

	void apply(std::complex<double>* inputState, std::complex<double>* outputState);
	void toMatrix(std::complex<double>** HMatrix);
	void print();
	std::string paramID();
};

#endif /* SZTRANSHAMILTONIAN_HPP_ */

