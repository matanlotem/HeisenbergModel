#include "szTransBasis.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "utils.hpp"
#include "szBasis.hpp"

// initialization
szTransBasis::szTransBasis(int n, int szUp, int k) {
	N = n;
	K = k;
	SzUp = szUp;
	szB = new szBasis(N, SzUp);

	combs = szB->getCombs();

	MapSzBasisIndex2TBasisIndex = new int[combs];
	MapSzBasisIndex2TNormFactor = new std::complex<double>[combs];
	for (int i=0; i<combs; i++) MapSzBasisIndex2TBasisIndex[i] = -1;

	if (SzUp == 0) {
		len = 1;
		MapSzBasisIndex2TBasisIndex[0] = 0;
		MapSzBasisIndex2TNormFactor[0] = 1;
		MapTBasisIndex2SzBasisState = new int[1];
		MapTBasisIndex2SzBasisState[0] = szB->getBasisState(0);
		MapTBasisIndex2TNormFactor = new double[1];
		MapTBasisIndex2TNormFactor[0] = 1;
	}
	else {
		len = 0;
		int tmpState;
		int tmpIndex;
		int* tmpSzBasisState = new int[N];
		std::complex<double>* tmpSzBasisStateValue = new std::complex<double>[N];
		int l;
		double norm;
		double eps = 0.000000000001; // for zero testing

		for (int i=0; i<combs; i++) {
			if (MapSzBasisIndex2TBasisIndex[i] == -1) {
				// for normalization factor
				for (int j=0; j<N; j++){
					tmpSzBasisState[j] = -1;
					tmpSzBasisStateValue[j] = 0;
				}

				tmpIndex = i;
				tmpState = szBasisIndex2szBasisState(tmpIndex);
				for (int j=0; j<N; j++) {
					// map szBasisIndex to TBasisIndex and TNormFactor (just phase here)
					MapSzBasisIndex2TBasisIndex[tmpIndex] = len;
					MapSzBasisIndex2TNormFactor[tmpIndex] = std::polar(1.0,2*M_PI * double(K * j) / N);

					// for normalization factor
					l = 0;
					while (tmpSzBasisState[l] != tmpState && tmpSzBasisState[l] != -1) l++;
					tmpSzBasisState[l] = tmpState;
					tmpSzBasisStateValue[l] += MapSzBasisIndex2TNormFactor[tmpIndex];

					// next state
					//tmpState = (tmpState >> 1) + ((tmpState & 1)<<(N-1)); // rotate right
					tmpState = ((tmpState << 1) & ~(1<<N)) + (tmpState >> (N-1)); // rotate left
					tmpIndex = szBasisState2szBasisIndex(tmpState);
				}

				// calculate normalization factor
				norm = 0;
				for (l=0; l<N; l++) norm += real(tmpSzBasisStateValue[l] * conj(tmpSzBasisStateValue[l]));
				norm = sqrt(norm);

				tmpIndex = i;
				tmpState = szBasisIndex2szBasisState(tmpIndex);

				do {
					if (norm < eps) // is zero => not a basis vector
						MapSzBasisIndex2TNormFactor[tmpIndex] = 0;
					else
						MapSzBasisIndex2TNormFactor[tmpIndex] = MapSzBasisIndex2TNormFactor[tmpIndex] / norm;
					tmpState = ((tmpState << 1) & ~(1<<N)) + (tmpState >> (N-1)); // rotate left
					tmpIndex = szBasisState2szBasisIndex(tmpState);
				//}
				} while (tmpIndex != i);

				if (norm >= eps) len++; // norm is not zero => is a basis vector
			}
		}

		delete[] tmpSzBasisState;
		delete[] tmpSzBasisStateValue;

		MapTBasisIndex2SzBasisState = new int[len];
		MapTBasisIndex2TNormFactor = new double[len];
		for (int i=0; i<len; i++) MapTBasisIndex2SzBasisState[i] = -1;
		for (int i=0; i<combs; i++) {
			if (MapTBasisIndex2SzBasisState[MapSzBasisIndex2TBasisIndex[i]] == -1 &&
					sqrt(real(MapSzBasisIndex2TNormFactor[i] * conj(MapSzBasisIndex2TNormFactor[i]))) >= eps) {
				MapTBasisIndex2SzBasisState[MapSzBasisIndex2TBasisIndex[i]] = szBasisIndex2szBasisState(i);
				MapTBasisIndex2TNormFactor[MapSzBasisIndex2TBasisIndex[i]] = abs(MapSzBasisIndex2TNormFactor[i]);
			}
		}
	}

}

szTransBasis::~szTransBasis() {
	delete[] MapSzBasisIndex2TBasisIndex;
	delete[] MapSzBasisIndex2TNormFactor;
	delete[] MapTBasisIndex2SzBasisState;
	delete[] MapTBasisIndex2TNormFactor;
	delete szB;
}


// getters
int szTransBasis::getCombs() {return combs;}
int szTransBasis::getLen() {return len;}
int szTransBasis::getN() {return N;}
int szTransBasis::getK() {return K;}
int szTransBasis::getSzUp() {return SzUp;}


int szTransBasis::szBasisIndex2TBasisIndex(int szBasisIndex) {
	if (szBasisIndex >=0 && szBasisIndex < combs)
		return MapSzBasisIndex2TBasisIndex[szBasisIndex];
	else
		return -1;
}

std::complex<double> szTransBasis::szBasisIndex2TNormFactor(int szBasisIndex) {
	if (szBasisIndex >=0 && szBasisIndex < combs)
		return MapSzBasisIndex2TNormFactor[szBasisIndex];
	else
		return 0;
}

int szTransBasis::TBasisIndex2szBasisState(int TBasisIndex) {
	if (TBasisIndex >= 0 && TBasisIndex < len)
		return MapTBasisIndex2SzBasisState[TBasisIndex];
	else
		return -1;
}

int szTransBasis::szBasisState2szBasisIndex(int szBasisState) {
	return szB->getBasisStateIndex(szBasisState);
}

int szTransBasis::szBasisIndex2szBasisState(int szBasisIndex) {
	return szB->getBasisState(szBasisIndex);
}

double szTransBasis::TBasisIndex2TNormFactor(int TBasisIndex) {
	return MapTBasisIndex2TNormFactor[TBasisIndex];
}

// states
void szTransBasis::newState(std::complex<double>* state) {
	for (int i = 0; i < len; i++) state[i] = 0;
}

void szTransBasis::newState(std::complex<double>* state, int TBasisIndex) {
	newState(state, TBasisIndex, std::complex<double>(1));
}

void szTransBasis::newState(std::complex<double>* state, int TBasisIndex, std::complex<double> TNormFactor) {
	newState(state);
	state[TBasisIndex] = TNormFactor;
}

void szTransBasis::copyState(std::complex<double>* inputState, std::complex<double>* outputState) {
	for (int i = 0; i < len; i++)
		outputState[i] = inputState[i];
}

std::complex<double> szTransBasis::scalarProd(std::complex<double>* state1, std::complex<double>* state2) {
	std::complex<double> prod = 0;
	for (int i = 0; i < len; i++)
		prod += state1[i] * conj(state2[i]);
	return prod;
}

double szTransBasis::getStateNorm(std::complex<double>* state) {
	return sqrt(real(scalarProd(state, state)));
}

double szTransBasis::normalizeState(std::complex<double>* state) {
	double norm = getStateNorm(state);
	for (int i = 0; i < len; i++)
		state[i] = state[i] / norm;
	return norm;
}

void szTransBasis::printState(std::complex<double>* state) {
	for (int i = 0; i < len; i++)
		printf("%+.3f%+.3fi\t", real(state[i]), imag(state[i]));
	printf("\n");
}
