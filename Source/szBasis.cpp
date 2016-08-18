#include "szBasis.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "utils.hpp"


// initialization
szBasis::szBasis(int n, int szUp) {
	N = n;
	SzUp = szUp;
	combs = utils::nCr(N,SzUp);
	statesArray = new int[combs];

	if (SzUp == 0)
		statesArray[0] = 0;
	else {
		// initializes stateNumbers
		int* stateNumbers = new int[SzUp];
		for (int i = 0; i < SzUp; i++) {
			stateNumbers[i] = i;
		}

		// insert first state
		int i = 0;
		statesArray[i++] = stateNumbers2State(stateNumbers);

		int j = SzUp - 1;
		while (stateNumbers[0] < (N - SzUp)) { // before last state is reached
			if (stateNumbers[j] < N - (SzUp - j)) {
				stateNumbers[j] ++;
				while (j + 1 < SzUp) {
					stateNumbers[j+1] = stateNumbers[j] + 1;
					j++;
				}
				statesArray[i++] = stateNumbers2State(stateNumbers);
			}
			else {
				j--;
			}
		}
		delete[] stateNumbers;

		// reverse list
		int tmp;
		for (i = 0; i < combs/2; i++) {
			tmp = statesArray[combs - 1 - i];
			statesArray[combs - 1 - i] = statesArray[i];
			statesArray[i] = tmp;
		}
	}

}

szBasis::~szBasis() {
	delete[] statesArray;
}

int szBasis::stateNumbers2State(int* stateNumbers) {
	int state = 0;
	for (int i = 0; i < SzUp; i++) {
		state += (1 << (N-stateNumbers[i]-1));
	}
	return state;
}

// getters
int szBasis::getCombs() {return combs;}
int szBasis::getLen() {return combs;}
int szBasis::getN() {return N;}
int szBasis::getSzUp() {return SzUp;}

int szBasis::getBasisState(int index) {
	return statesArray[index];
}

int szBasis::getBasisStateIndex(int basisState) {
	// binary search for state in state array. if not found return -1
	int l = 0;
	int r = combs - 1;
	int m;
	int i = 0;
	while (l < r) {
		m = (l + r) / 2;
		if (basisState == statesArray[m]) {
			r = m;
			l = m;
		}
		else if (basisState < statesArray[m]) {
			r = m - 1;
		}
		else {
			l = m + 1;
		}
		i++;
	}

	if (statesArray[l] == basisState)
		return l;
	else
		return -1;
}

// states
void szBasis::newState(double* state) {
	for (int i = 0; i < combs; i++)
		state[i] = 0;
}

void szBasis::newState(double* state, int basisIndex) {
	newState(state);
	state[basisIndex] = 1;
}

void szBasis::copyState(double* inputState, double* outputState) {
	for (int i = 0; i < getLen(); i++)
		outputState[i] = inputState[i];
}

double szBasis::scalarProd(double* state1, double* state2) {
	double prod = 0;
	for (int i = 0; i < getLen(); i++)
		prod += state1[i] * state2[i];
	return prod;
}

double szBasis::getStateNorm(double* state) {
	return sqrt(scalarProd(state, state));
}

double szBasis::normalizeState(double* state) {
	double norm = getStateNorm(state);
	for (int i = 0; i < getLen(); i++)
		state[i] = state[i] / norm;
	return norm;
}

void szBasis::printState(double* state) {
	for (int i = 0; i < getLen(); i++)
		printf("%+.3f\t", state[i]);
	printf("\n");
}
