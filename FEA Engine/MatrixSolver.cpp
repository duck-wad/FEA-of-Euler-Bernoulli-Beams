#include "MatrixSolver.h"
#include "Utils.h"

std::vector<double> GaussSeidel(const std::vector<std::vector<double>>& A, const std::vector<double>& b, const double errorTol, const int maxIter) {
	if (A.size() != b.size() || A[0].size() != b.size()) {
		throw std::invalid_argument("Matrix and vector must have same dimensions");
	}
	//initialize guess to zero
	std::vector<double> x(A.size(), 0.0);
	int numIter = 0;

	for (size_t iter = 0; iter < maxIter; iter++) {
		numIter = iter + 1;
		//store results from previous iteration
		std::vector<double> x_old = x;
		//compute the solution for each row sequentially starting with the first row
		for (size_t i = 0; i < A.size(); i++) {
			double sum = 0.0;
			for (size_t j = 0; j < i; j++) {
				sum += A[i][j] * x[j];
			}
			for (size_t j = i + 1; j < A.size(); j++) {
				sum += A[i][j] * x_old[j];
			}
			//update result for current iteration
			x[i] = (b[i] - sum) / A[i][i];
		}
		//check convergence. error for each value must be below tolerance
		bool errorFlag = false;
		for (size_t i = 0; i < A.size(); i++) {
			double error = (x[i] - x_old[i]) / x[i];
			if (error > errorTol) errorFlag = true;
		}
		if (!errorFlag) break;

	}

	if (numIter == maxIter) {
		throw std::runtime_error("Gauss-Seidel method failed to converge within maximum iterations");
	}

	return x;
}

