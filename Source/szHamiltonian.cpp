/*
 * Classes.cpp
 *
 *  Created on: Jul 30, 2016
 *      Author: Matan
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "utils.hpp"
#include "szHamiltonian.hpp"
#include "szBasis.hpp"

szHamiltonian::szHamiltonian(double jxy, double jz, double hz, bool cyc, szBasis * b):
	Jxy(jxy), Jz(jz), Hz(hz), cyclic(cyc) ,basis(b) {

	N = basis->getN();
	SzUp = basis->getSzUp();
	len = basis->getLen();

	stateValues = new double[len];
	ladderMatrix = new int*[len];
	ladderMatrix[0] = new int[len*2*N];
	for (int i = 1; i < len; i++) ladderMatrix[i] = ladderMatrix[i - 1] + 2*N;

	// Precalc Hamiltonian
	double stateValue;
	int currentBasisState;

	for (int i = 0; i < len; i++) {
		stateValue = 0;
		currentBasisState = basis->getBasisState(i);
		for (int j = 1; j <= N - (!cyclic); j++) { // if cyclic do N iterations, otherwise do N-1 iterations
			ladderMatrix[i][j-1] = basis->getBasisStateIndex(ladderUpDown(j, currentBasisState));
			ladderMatrix[i][N+j-1] = basis->getBasisStateIndex(ladderDownUp(j, currentBasisState));

			stateValue += (Jz * spinValue(j, currentBasisState) * spinValue(j % N + 1, currentBasisState));
			stateValue += Hz * spinValue(j, currentBasisState);
		}
		if (!cyclic) // if not cyclic get last spin
			stateValue += Hz * spinValue(N, currentBasisState);
		stateValues[i] = stateValue;
	}

}

szHamiltonian::~szHamiltonian() {
	delete[] stateValues;
	delete[] ladderMatrix[0];
	delete[] ladderMatrix;
}


int szHamiltonian::ladderUp(int i, int basisState) {
	if (!((basisState >> (N-i)) & 1)) // if bit i is 0
		return  (basisState ^ (1 << (N - i))); // flip it
	return -1;
}

int szHamiltonian::ladderDown(int i, int basisState) {
	if ((basisState >> (N - i)) & 1) // if bit i is 1
		return  (basisState ^ (1 << (N - i))); // flip it
	return -1;
}

int szHamiltonian::ladderUpDown(int i, int basisState) {
	int newState = ladderUp(i, basisState);
	if (newState != -1)
		newState = ladderDown(i % N + 1, newState);
	return newState;
}

int szHamiltonian::ladderDownUp(int i, int basisState) {
	int newState = ladderDown(i, basisState);
	if (newState != -1)
		newState = ladderUp(i % N + 1, newState);
	return newState;
}

double szHamiltonian::spinValue(int i, int basisState) {
	return ((basisState >> (N-i)) & 1) - 0.5;
}

void szHamiltonian::apply(double* inputState, double* outputState) {
	basis->newState(outputState);
	int outputBasisStateIndex;

	for (int i = 0; i < len; i++) {
		if (inputState[i] != 0) {
			for (int j = 1; j <= N - (!cyclic); j++) { // if cyclic do N iterations, otherwise do N-1 iterations
				outputBasisStateIndex = ladderMatrix[i][j - 1];
				if (outputBasisStateIndex >= 0 && outputBasisStateIndex < len)
					outputState[outputBasisStateIndex] += inputState[i] * Jxy / 2;
				outputBasisStateIndex = ladderMatrix[i][N + j - 1];
				if (outputBasisStateIndex >= 0 && outputBasisStateIndex < len)
					outputState[outputBasisStateIndex] += inputState[i] * Jxy / 2;
			}
			outputState[i] += inputState[i] * stateValues[i];
		}

	}
}

void szHamiltonian::toMatrix(double** HMatrix) {
	double * tmpState = new double[basis->getLen()];
	for (int i = 0; i < len; i++) {
		basis->newState(tmpState,i);
		apply(tmpState, HMatrix[i]);
	}
	delete[] tmpState;
}

void szHamiltonian::print() {
	double** HMatrix = new double*[len];
	HMatrix[0] = new double[len*len];
	for (int i = 1; i < len; i++) HMatrix[i] = HMatrix[i - 1] + len;

	toMatrix(HMatrix);
	utils::printMatrix(len,HMatrix);

	delete[] HMatrix[0];
	delete[] HMatrix;
}

std::string szHamiltonian::paramID() {
	char buffer [50];
	sprintf(buffer,"N%d_SzUp%d_Jxy%.1f_Jz%.1f_Hz%.1f_cyc%d", basis->getN(), basis->getSzUp(), Jxy, Jz, Hz, cyclic);
	std::string paramStr(buffer);
	return paramStr;
}
