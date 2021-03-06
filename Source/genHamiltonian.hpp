/*
 * genHamiltonian.hpp
 *
 *  Created on: Aug 17, 2016
 *      Author: Matan
 */

#ifndef SOURCE_GENHAMILTONIAN_HPP_
#define SOURCE_GENHAMILTONIAN_HPP_

template <class T>
class genHamiltonian {
public:
	virtual ~genHamiltonian();
	virtual void apply(T* inputState, T* outputState) = 0;
	virtual void toMatrix(T** HMatrix) = 0;
	virtual void print() = 0;
};

template <class T> genHamiltonian<T>::~genHamiltonian() {}

#endif /* SOURCE_GENHAMILTONIAN_HPP_ */
