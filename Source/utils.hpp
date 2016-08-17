/*
 * utils.hpp
 *
 *  Created on: Jul 30, 2016
 *      Author: Matan
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_
#include <string>
#include <complex>

class utils {
public:
	static long nCr(int n, int r);
	static void printBinary(int value);
	static void printBinary(int value, int len);

	static void printMatrix(int dim, int** matrix);
	static void printMatrix(int rows, int cols, int** matrix);
	static void printMatrix(int dim, double** matrix);
	static void printMatrix(int rows, int cols, double** matrix);
	static void printMatrix(int dim, std::complex<double>** matrix);
	static void printMatrix(int rows, int cols, std::complex<double>** matrix);

	static void copyMatrix(int dim, double** inputMatrix, double** outputMatrix);
	static void copyMatrix(int dim, std::complex<double>** inputMatrix, std::complex<double>** outputMatrix);

	static int diagnolize(int dim, double **matrix, double *EigenValues);
	static int diagnolize(int dim, double** matrix, double* EigenValues, std::string action);
	static int diagnolize(int dim, std::complex<double>** matrix, double* EigenValues);
	static int diagnolize(int dim, std::complex<double>** matrix, double* EigenValues, std::string action);

	static void clearLog(const char * logfile);
	static void logprintf(const char * logfile, const char * format, ...);
};


#endif /* UTILS_HPP_ */
