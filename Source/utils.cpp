/*
 * utils.cpp
 *
 *  Created on: Jul 30, 2016
 *      Author: Matan
 */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <complex>
#include <cstdarg>

#include "utils.hpp"

extern "C" void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda,
	double* w, double* work, int* lwork, int* info);

extern "C" void zheev_(char* jobz, char* uplo, int* n, std::complex<double>* a, int* lda,
	double* w, std::complex<double>* work, int* lwork, double* rwork, int* info);

long utils::nCr(int n, int r) {
	if (r == 0)
		return 1;

	if (n - r < r)
		return nCr(n, n - r);

	long long res = 1;

	int* rs = new int[r];
	for (int i = 0; i < r; i++)
		rs[i] = 1;

	for (int i = n - r + 1; i <= n; i++) {
		res *= i;
		for (int j = 1; j <= r; j++) {
			if (rs[j - 1] && res % j == 0) {
				res /= j;
				rs[j - 1] = 0;
			}
		}
	}

	delete[] rs;
	return res;
}

void utils::printBinary(int value) {
	printBinary(value, 32);
}

void utils::printBinary(int value, int len) {
	for (int i=len-1; i>=0; i--)
		printf("%d",(value>>i & 1));
	printf("\n");
}

void utils::printMatrix(int dim, int** matrix) {
	printMatrix(dim, dim, matrix);
}

void utils::printMatrix(int rows, int cols, int** matrix) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) printf("%d\t",matrix[i][j]);
		printf("\n\n");
	}
}

void utils::printMatrix(int dim, double** matrix) {
	printMatrix(dim, dim, matrix);
}

void utils::printMatrix(int rows, int cols, double** matrix) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) printf("%+.2f\t",matrix[i][j]);
		printf("\n\n");
	}
}

void utils::printMatrix(int dim, std::complex<double>** matrix) {
	printMatrix(dim, dim, matrix);
}

void utils::printMatrix(int rows, int cols, std::complex<double>** matrix) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) printf("%+.2f%+.2fi\t",real(matrix[i][j]),imag(matrix[i][j]));
		printf("\n\n");
	}
}

void utils::copyMatrix(int dim, double** inputMatrix, double** outputMatrix) {
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++) outputMatrix[i][j] = inputMatrix[i][j];
}

void utils::copyMatrix(int dim, std::complex<double>** inputMatrix, std::complex<double>** outputMatrix) {
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++) outputMatrix[i][j] = inputMatrix[i][j];
}


int utils::diagnolize(int dim, double **matrix, double *EigenValues) {
	return diagnolize(dim, matrix, EigenValues, "N");
}

int utils::diagnolize(int dim, double **matrix, double *EigenValues, std::string action) {
	double* workMatrix;
	char* cmnd = "N";
	if (action == "V")
		cmnd = "V";

	if (action == "C") {
		workMatrix = new double[dim*dim];
		for (int i=0; i < dim; i++)
			for (int j=0; j<dim; j++) workMatrix[i*dim+j] = matrix[i][j];
	}
	else
		workMatrix = matrix[0];

	int n = dim, lda = dim;
	int info = 0, lwork;
	double wkopt;

	lwork = -1;
	dsyev_(cmnd, "U", &n, workMatrix, &lda, EigenValues, &wkopt, &lwork, &info);
	lwork = (int)wkopt;
	double *work = new double[lwork];
	dsyev_(cmnd, "U", &n, workMatrix, &lda, EigenValues, work, &lwork, &info);

	return info;
}


int utils::diagnolize(int dim, std::complex<double>** matrix, double* EigenValues) {
	return diagnolize(dim, matrix, EigenValues, "N");
}

int utils::diagnolize(int dim, std::complex<double>** matrix, double* EigenValues, std::string action) {
	std::complex<double>* workMatrix;
	char* cmnd = "N";
	if (action == "V")
		cmnd = "V";

	if (action == "C") {
		workMatrix = new std::complex<double>[dim*dim];
		for (int i=0; i < dim; i++)
			for (int j=0; j<dim; j++) workMatrix[i*dim+j] = matrix[i][j];
	}
	else
		workMatrix = matrix[0];

	int n = dim, lda = dim;
	int info = 0, lwork;
	std::complex<double> wkopt;
	double* rwork = new double[std::max(1,3*n-2)];

	lwork = -1;
	zheev_(cmnd, "U", &n, workMatrix, &lda, EigenValues, &wkopt, &lwork, rwork, &info);
	lwork = (int)real(wkopt);
	std::complex<double>* work = new std::complex<double>[lwork];
	zheev_(cmnd, "U", &n, workMatrix, &lda, EigenValues, work, &lwork, rwork, &info);

	return info;
}

void utils::clearLog(const char * logfile) {
	FILE * file = fopen(logfile, "w");
	fprintf (file, "");
	fclose(file);
}

void utils::logprintf(const char * logfile, const char * format, ...) {
	FILE * file = fopen(logfile, "a");
	va_list args;
	va_start (args, format);
	vfprintf (file, format, args);
	va_end (args);
	fclose(file);
}
