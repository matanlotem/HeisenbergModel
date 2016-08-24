/*
 * genBasis.hpp
 *
 *  Created on: Aug 17, 2016
 *      Author: Matan
 */

#ifndef SOURCE_GENBASIS_HPP_
#define SOURCE_GENBASIS_HPP_


template <class T>
class genBasis {
public:
	virtual ~genBasis();
	virtual int getLen() = 0;
	virtual void newState(T* state) = 0;
	virtual void newState(T* state, int ) = 0;
	virtual void copyState(T* inputState, T* outputState) = 0;
	virtual T scalarProd(T* state1, T* state2) = 0;
	virtual double getStateNorm(T* state) = 0;
	virtual void printState(T* state) = 0;
};

template <class T> genBasis<T>::~genBasis() {}

#endif /* SOURCE_GENBASIS_HPP_ */
