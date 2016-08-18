/*
 * szBasis.hpp
 *
 *  Created on: Jul 30, 2016
 *      Author: Matan
 */

#ifndef SZBASIS_HPP_
#define SZBASIS_HPP_
#include "genBasis.hpp"

class szBasis : public genBasis<double>{
private:
	int N;
	int SzUp;
	int combs;
	int* statesArray;

	int stateNumbers2State(int* stateNumbers);
public:
	szBasis(int n, int szUp);
	~szBasis();
	int getCombs();
	int getLen();
	int getN();
	int getSzUp();
	int getBasisState(int index);
	int getBasisStateIndex(int basisState);

	void newState(double* state);
	void newState(double* state, int basisIndex);
	void copyState(double* inputState, double* outputState);
	double scalarProd(double* state1, double* state2);
	double getStateNorm(double* state);
	double normalizeState(double* state);
	void printState(double* state);
};



#endif /* SZBASIS_HPP_ */
