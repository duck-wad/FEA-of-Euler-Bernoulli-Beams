#include "MatrixSolver.h"
#include "Utils.h"

std::vector<double> GaussSeidel(const std::vector<std::vector<double>>& A, const std::vector<double>& b, const double errorTol, const int maxIter) {
	if (A.size() != b.size() || A[0].size() != b.size()) {
		throw std::invalid_argument("Matrix and vector must have same dimensions");
	}
	//initialize guess to zero
	std::vector<double> x(A.size(), 0.0);
	int numIter = 0;

	for (int iter = 0; iter < maxIter; iter++) {
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

std::vector<double> GaussianElimination(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {

	if (A.size() != b.size() || A[0].size() != b.size()) {
		throw std::invalid_argument("Matrix and vector must have same dimensions");
	}

	std::vector<std::vector<double>> A_REF = A;
	std::vector<double> b_REF = b;

	ForwardElimination(A_REF, b_REF);

	//check if A_REF is singular (zero exists on the diagonal)
	if (isSingular(A_REF)) {
		throw std::invalid_argument("Matrix is singular, system cannot be solved for a unique solution");
	}

	std::vector<double> x = BackwardSubstitution(A_REF, b_REF);
	return x;
}

//if L needs to be filled, the function expects L to be sized before being passed in
void ForwardElimination(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<std::vector<double>>& L) {

	if (A.size() != b.size() || A[0].size() != b.size()) {
		throw std::invalid_argument("Matrix and vector must have same dimensions");
	}

	bool fillL = !L.empty();
	if (fillL) {
		if (L.size() != A.size() || L[0].size() != A[0].size()) {
			throw std::invalid_argument("L is not sized correctly");
		}
		for (size_t i = 0; i < L.size(); i++) {
			L[i][i] = 1.0;
		}
	}

	//append b to last column of A matrix to put into augmented form to streamline row operations
	std::vector<std::vector<double>> temp = transpose(A);
	temp.push_back(b);
	std::vector<std::vector<double>> aug_A = transpose(temp);

	//perform row operations to bring aug_A to REF
	for (size_t piv = 0; piv < aug_A.size(); piv++) {
		for (size_t i = piv + 1; i < aug_A.size(); i++) {
			double f = aug_A[i][piv] / aug_A[piv][piv];
			std::vector<double> pivRow = aug_A[piv] * f;
			aug_A[i] = aug_A[i] - pivRow;
			if (fillL) {
				L[i][piv] = f;
			}
		}
	}
	//set b to last column of augmented matrix
	//un-augment to get the REF for A
	temp = transpose(aug_A);
	b = temp[A.size()];
	temp.pop_back();
	A = transpose(temp);
}

std::vector<double> BackwardSubstitution(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {

	std::vector<double> x(b.size());
	for (size_t i = A.size(); i-- > 0; ) {
		if (i == A.size() - 1) {
			x[i] = b[i] / A[i][i];
		}
		else {
			double sum = 0.0;
			for (size_t j = i + 1; j < A.size(); j++) {
				sum += A[i][j] * x[j];
			}
			x[i] = (b[i] - sum) / A[i][i];
		}
	}
	return x;
}

//if L is not passed in, dummy function is called which then calls actual function with an empty L
void ForwardElimination(std::vector<std::vector<double>>& A, std::vector<double>& b) {
	std::vector<std::vector<double>> L = {};
	ForwardElimination(A, b, L);
}






