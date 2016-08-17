/*
 * genBasis.hpp
 *
 *  Created on: Aug 17, 2016
 *      Author: Matan
 */

#ifndef SOURCE_GENBASIS_HPP_
#define SOURCE_GENBASIS_HPP_


class genBasis {
public:
	int getLen();
	void copyMixedState(double* inputState, double* outputState);
	double scalarProd(double* state1, double* state2);
	double getMixedStateNorm(double* state);
};


#endif /* SOURCE_GENBASIS_HPP_ */
