/*
 * szHamiltonian.hpp
 *
 *  Created on: Jul 30, 2016
 *      Author: Matan
 */

#ifndef SZHAMILTONIAN_HPP_
#define SZHAMILTONIAN_HPP_
#include <string>
#include "genHamiltonian.hpp"
#include "szBasis.hpp"

class szHamiltonian : public genHamiltonian<double> {
protected:
	int N;
	int SzUp;
	int len;
	double Jxy;
	double Jz;
	double Hz;
	bool cyclic;
	szBasis * basis;
private:
	int** ladderMatrix;
	double* stateValues;
public:
	szHamiltonian(double jxy, double jz, double hz, bool cyc, szBasis * b);
	~szHamiltonian();

	int ladderUp(int i, int basisState);
	int ladderDown(int i, int basisState);
	int ladderUpDown(int i, int basisState);
	int ladderDownUp(int i, int basisState);
	double spinValue(int i, int basisState);

	void apply(double* inputState, double* outputState);
	void toMatrix(double** HMatrix);
	void print();
	std::string paramID();
};




#endif /* SZHAMILTONIAN_HPP_ */
