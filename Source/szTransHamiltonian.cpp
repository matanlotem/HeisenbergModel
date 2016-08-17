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
#include "szTransHamiltonian.hpp"
#include "szTransBasis.hpp"

szTransHamiltonian::szTransHamiltonian(double jxy, double jz, double hz, szTransBasis * b):
	Jxy(jxy), Jz(jz), Hz(hz) ,basis(b) {

	N = basis->getN();
	K = basis->getK();
	SzUp = basis->getSzUp();
	len = basis->getLen();

	stateValues = new std::complex<double>[len];
	ladderMatrix = new int*[len];
	ladderMatrix[0] = new int[len*2*N];
	for (int i = 1; i < len; i++) ladderMatrix[i] = ladderMatrix[i - 1] + 2*N;

	ladderValue = new std::complex<double>*[len];
	ladderValue[0] = new std::complex<double>[len*2*N];
	for (int i = 1; i < len; i++) ladderValue[i] = ladderValue[i - 1] + 2*N;

	// Precalc Hamiltonian
	std::complex<double> stateValue;
	int currentSzBasisState;
	int szBasisIndex;
	double norm;

	for (int i = 0; i < len; i++) {
		stateValue = 0;
		currentSzBasisState = basis->TBasisIndex2szBasisState(i);
		for (int j = 1; j <= N; j++) {
			szBasisIndex = basis->szBasisState2szBasisIndex(szLadderUpDown(j,currentSzBasisState));
			ladderMatrix[i][j-1] = basis->szBasisIndex2TBasisIndex(szBasisIndex);
			ladderValue [i][j-1] = 0;
			if (ladderMatrix[i][j-1] != -1 && basis->szBasisIndex2TNormFactor(szBasisIndex) != 0.0) { // zero division protection
				ladderValue [i][j-1] = basis->TBasisIndex2TNormFactor(i) / basis->szBasisIndex2TNormFactor(szBasisIndex) * (Jxy / 2.0);
			}

			szBasisIndex = basis->szBasisState2szBasisIndex(szLadderDownUp(j,currentSzBasisState));
			ladderMatrix[i][N+j-1] = basis->szBasisIndex2TBasisIndex(szBasisIndex);
			ladderValue [i][N+j-1] = 0;
			if (ladderMatrix[i][N+j-1] != -1 && basis->szBasisIndex2TNormFactor(szBasisIndex) != 0.0) { // zero division protection
				ladderValue [i][N+j-1] = basis->TBasisIndex2TNormFactor(i) / basis->szBasisIndex2TNormFactor(szBasisIndex) * (Jxy / 2.0);
			}

			stateValue += (Jz * szSpinValue(j, currentSzBasisState) * szSpinValue(j % N + 1, currentSzBasisState));
			stateValue += Hz * szSpinValue(j, currentSzBasisState);
		}

		stateValues[i] = stateValue;
	}

	/*printf("ladderMatrix:\n");
	utils::printMatrix(len, 2*N, ladderMatrix);
	printf("\nladderValue:\n");
	utils::printMatrix(len, 2*N, ladderValue);/**/
}

szTransHamiltonian::~szTransHamiltonian() {
	delete[] stateValues;
	delete[] ladderMatrix[0];
	delete[] ladderMatrix;
	delete[] ladderValue[0];
	delete[] ladderValue;
}

// szBasis operators (should be moved to szBasisStates class)
int szTransHamiltonian::szLadderUp(int i, int basisState) {
	if (!((basisState >> (N-i)) & 1)) // if bit i is 0
		return  (basisState ^ (1 << (N - i))); // flip it
	return -1;
}

int szTransHamiltonian::szLadderDown(int i, int basisState) {
	if ((basisState >> (N - i)) & 1) // if bit i is 1
		return  (basisState ^ (1 << (N - i))); // flip it
	return -1;
}

int szTransHamiltonian::szLadderUpDown(int i, int basisState) {
	int newState = szLadderUp(i, basisState);
	if (newState != -1)
		newState = szLadderDown(i % N + 1, newState);
	return newState;
}

int szTransHamiltonian::szLadderDownUp(int i, int basisState) {
	int newState = szLadderDown(i, basisState);
	if (newState != -1)
		newState = szLadderUp(i % N + 1, newState);
	return newState;
}

double szTransHamiltonian::szSpinValue(int i, int basisState) {
	return ((basisState >> (N-i)) & 1) - 0.5;
}

void szTransHamiltonian::apply(std::complex<double>* inputState, std::complex<double>* outputState) {
	basis->newState(outputState);
	int outputTBasisIndex;

	for (int i = 0; i < len; i++) {
		if (inputState[i] != std::complex<double>(0)) {
			for (int j = 1; j <= N; j++) {
				outputTBasisIndex = ladderMatrix[i][j - 1];

				if (outputTBasisIndex >= 0 && outputTBasisIndex < len)
					outputState[outputTBasisIndex] += ladderValue[i][j - 1] * inputState[i];

				outputTBasisIndex = ladderMatrix[i][N + j - 1];
				if (outputTBasisIndex >= 0 && outputTBasisIndex < len)
					outputState[outputTBasisIndex] += ladderValue[i][N + j - 1] * inputState[i];
			}

			outputState[i] += inputState[i] * stateValues[i];
		}

	}
}

void szTransHamiltonian::toMatrix(std::complex<double>** HMatrix) {
	std::complex<double> * tmpState = new std::complex<double>[basis->getLen()];
	for (int i = 0; i < len; i++) {
		basis->newState(tmpState,i);
		apply(tmpState, HMatrix[i]);
	}
	delete[] tmpState;
}

void szTransHamiltonian::print() {
	std::complex<double>** HMatrix = new std::complex<double>*[len];
	HMatrix[0] = new std::complex<double>[len*len];
	for (int i = 1; i < len; i++) HMatrix[i] = HMatrix[i - 1] + len;

	toMatrix(HMatrix);
	utils::printMatrix(len,HMatrix);

	delete[] HMatrix[0];
	delete[] HMatrix;
}

std::string szTransHamiltonian::paramID() {
	char buffer [50];
	sprintf(buffer,"N%d_K%d_SzUp%d_Jxy%.1f_Jz%.1f_Hz%.1f", basis->getN(), basis->getK(), basis->getSzUp(), Jxy, Jz, Hz);
	std::string paramStr(buffer);
	return paramStr;
}
