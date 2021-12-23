/* FILE NAME   : main.cpp
 * PROGRAMMER  : Kononov Svyatoslav '4117.
 * LAST UPDATE : 26.11.2021.
 */

#include <iostream>
#include <vector>
#include <set>
#include <complex>
#include "Eigen/Eigenvalues"

// Use linear algebra library namespace
using namespace Eigen;

// Epsilon constant
const double Eps = 1e-13;

/* Compare number with epsilon function.
 * ARGUMENTS:
 *     - number:
 *         std::complex<double> c;
 * RETURNS:
 *     (bool) true if number modulus less than epsilon, false otherwise.
 */
bool isZero(std::complex<double> c) {
	return c.real() < Eps && c.real() > -Eps;
} /* End of 'matRound' function */

/* Round matrix elements function.
 * ARGUMENTS:
 *     - matrix:
 *         const MatrixXcd &M;
 * RETURNS:
 *     (MatrixXcd) rounded matrix.
 */
MatrixXcd matRound(MatrixXcd &M) {
	// Create rounded matrix
	MatrixXcd cpy = M;

	// Round all elements
	for (uint32_t i = 0; i < cpy.rows(); i++)
		for (uint32_t j = 0; j < cpy.cols(); j++)
			if (isZero(cpy(i, j)))
				cpy(i, j) = 0;

	// Return rounded matrix
	return cpy;
} /* End of 'matRound' function */

/* Get submatrix function.
 * ARGUMENTS:
 *     - matrix:
 *         const MatrixXd &M;
 *     - rows/columns submatrix set:
 *         const std::set<uint32_t> &sub;
 * RETURNS:
 *     (MatrixXd) submatrix.
 */
MatrixXd submatrix(const MatrixXd &M, const std::set<uint32_t> &sub = {}) {
	// Create submatrix
	MatrixXd res(sub.size(), sub.size());

	// Copy matrix elements of submatrix rows/columns set
	uint32_t i = 0, j = 0;
	for (auto x : sub) {
		for (auto y : sub) {
			res.block(i, j, 1, 1) = M.block(x, y, 1, 1);
			j++;
		}
		i++;
		j = 0;
	}

	// Return submatrix
	return res;
} /* End of 'submatrix' function */

/* Calculate matrix invariant function.
 * ARGUMENTS:
 *     - matrix:
 *         const MatrixXd &M;
 *     - submatrix order:
 *         const uint32_t &k;
 *     - cycle starting position:
 *         const uint32_t &s;
 *     - submatrix rows and columns:
 *         const std::set<uint32_t> &sub;
 * RETURNS:
 *     (double) invariant.
 */
double invariant(const MatrixXd &M, const uint32_t &k, const uint32_t &s = 0, const std::set<uint32_t> &sub = {}) {
	// Check if matrix is empty
	assert((M.rows() || M.cols()) && "Empty matrix");

	// Check if matrix is square
	assert((M.rows() == M.cols()) && "Matrix is not square");

	// Check for submatrix order and return its determinant
	if (k == 0)
		// Return submatrix determinant
		return submatrix(M, sub).determinant();

	// Calculate sum of main diagonal minors
	double val = 0;
	for (uint32_t i = s; i < M.rows() - k + 1; i++) {
		std::set<uint32_t> subnew = sub;
		subnew.insert(i);
		val += invariant(M, k - 1, i + 1, subnew);
	}

	// Return invariant
	return val;
} /* End of 'invariant' function */

/* Calculate matrix invariants function.
 * ARGUMENTS:
 *     - matrix:
 *         const MatrixXd &M;
 * RETURNS:
 *     (std::vector<double>) invariants.
 */
std::vector<double> invariants(const MatrixXd &M) {
	// Create invariants vector and allocate memory for it
	std::vector<double> invs;
	invs.reserve(M.rows());

	// Calculate invariants for all submatrix orders
	for (uint32_t k = 1; k <= M.rows(); k++)
		invs.push_back(invariant(M, k));

	// Return invariants vector
	return invs;
} /* End of 'invariants' function */

/* Get equation from equation function.
 * ARGUMENTS:
 *     - matrix:
 *         const MatrixXd &N;
 * RETURNS: None.
 */
void outCurve(const MatrixXcd &N) {
	std::cout
		<< N(0, 0).real() << " x^2 + " << N(1, 1).real() << " y^2 + "
		<< N(0, 1).real() * 2. << " xy + "
		<< N(0, 2).real() * 2. << " x + " << N(1, 2).real() * 2. << " y + "
		<< N(2, 2).real() << " = 0" << std::endl;
} /* End of 'outCurves' function */

/* Get equation from equation function.
 * ARGUMENTS:
 *     - matrix:
 *         const MatrixXd &N;
 * RETURNS: None.
 */
void outSurface(const MatrixXcd &N) {
	std::cout
		<< N(0, 0).real() << " x^2 + " << N(1, 1).real() << " y^2 + " << N(2, 2).real() << " z^2 + "
		<< N(0, 1).real() * 2. << " xy + " << N(0, 2).real() * 2. << " xz + " << N(1, 2).real() * 2. << " yz + "
		<< N(0, 3).real() * 2. << " x + " << N(1, 3).real() * 2. << " y + " << N(2, 3).real() * 2. << " z + "
		<< N(3, 3).real() << " = 0" << std::endl;
} /* End of 'outSurface' function */

/* Determine 2nd order curve or surface function.
 * ARGUMENTS: None.
 * RETURNS: None.
 */
void CurvesAndSurfaces2ndOrder(void) {
	// Run flag
	bool flag = true;

	// Show equation flag
	bool showeq = false;

	// 2nd order coefficients
	double A, B, C, D, E, F, G, H, I, J;

	// 2nd order invariants
	std::vector<double> invs, subinvs;

	// UI variable
	char str;

	// Coefficients matrix
	MatrixXd coeffs;

	// Eigen vectors matrix
	MatrixXcd V;

	// New basis matrix
	MatrixXcd N;

	// Eigen solver
	EigenSolver<MatrixXd> es;

	// Set output parameters
	std::cout << std::scientific << std::fixed;

	// Run
	while (flag) {
		// Introduce possible variants to user
		std::cout << "--------------------------------\n";
		std::cout << "Choose what you want to explore:\n";
		std::cout << "0 - Quit\n";
		std::cout << "1 - Second order curves\n";
		std::cout << "2 - Second order surfaces\n";
		std::cout << "--------------------------------\n";

		// Get response
		std::cin >> str;

		// Run one of cases determined by user's response
		switch (str) {

		// 2nd order curve exploration variant
		case '1':

			// Introduce to user how to input 2nd order curve coefficients
			std::cout << "Enter the coefficients of the 2nd order curve as in the example:\n";
			std::cout << "If you have Ax^2 + By^2 + Cxy + Dx + Ey + F = 0\n";
			std::cout << "You should input <<A B C D E F>>\n";

			// Get 2nd order curve coefficients
			std::cin >> A >> B >> C >> D >> E >> F;

			// Divide repeating coefficients by 2
			C /= 2;
			D /= 2;
			E /= 2;

			// Construct coefficients matrix
			coeffs = MatrixXd(3, 3);
			coeffs <<
				A, C, D,
				C, B, E,
				D, E, F;

			// Calculate invariants
			subinvs = invariants(coeffs.block(0, 0, 2, 2));
			invs = invariants(coeffs);

			// Determine 2nd order curve by invariants combinations
			std::cout << std::endl;
			if (subinvs[1] != 0) {
				if (subinvs[1] > 0) {
					if (subinvs[0] * invs[2] < 0)
						std::cout << "Ellipse\n";
					else if (invs[2] == 0)
						std::cout << "Point\n";
					else if (subinvs[0] * invs[2] > 0)
						std::cout << "Imaginary ellipse\n";
				} else if (subinvs[1] < 0) {
					if (invs[2] != 0)
						std::cout << "Hyperbola\n";
					else
						std::cout << "Intersecting lines\n";
				}
			} else {
				if (invs[2] != 0)
					std::cout << "Parabola\n";
				else {
					if (invs[1] < 0)
						std::cout << "Parallel lines\n";
					else if (invs[1] == 0)
						std::cout << "Line\n";
					else
						std::cout << "Imaginary parallel lines\n";
				}
			}

			// Calculate eigen vectors of quadratic coefficients submatrix
			es = EigenSolver<MatrixXd>(coeffs.block(0, 0, 2, 2));

			// Prepare eigen vectors matrix
			V = Matrix3cd::Identity();
			V.block(0, 0, 2, 2) = es.eigenvectors();

			// Apply new basis
			N = V.inverse() * coeffs * V;

			// Round near zero elements
			N = matRound(N);

			// Show equation
			std::cout << std::endl;
			std::cout << "Applying new basis" << std::endl;
			outCurve(N);
			std::cout << std::endl;

			// Completing the square
			for (uint32_t i = 0; i < 2; i++)
				if (N(i, i) != std::complex<double>(0, 0) &&
					N(i, 2) != std::complex<double>(0, 0)) {
					N(2, 2) -= N(i, 2) * N(i, 2) / N(i, i);
					N(i, 2) = N(2, i) = 0;
					showeq = true;
				}

			// Show equation
			if (showeq) {
				std::cout << std::endl;
				std::cout << "Completing squares" << std::endl;
				outCurve(N);
				std::cout << std::endl;
			}

			// Reset flag
			showeq = false;

			// Orthogonal substitution
			if (N(0, 2) != std::complex<double>(0, 0) &&
				N(1, 2) != std::complex<double>(0, 0)) {
				std::complex<double> u = sqrt(
					N(0, 2) * N(0, 2) +
					N(1, 2) * N(1, 2));
				N(0, 2) = N(2, 0) = u;
				N(1, 2) = N(2, 1) = 0;
				showeq = true;
			}

			// Show equation
			if (showeq) {
				std::cout << std::endl;
				std::cout << "Applying orthogonal substitution" << std::endl;
				outCurve(N);
				std::cout << std::endl;
			}

			// Reset flag
			showeq = false;

			// Nullify free temp if at least one linear coefficient is not a zero
			if (N(0, 2) != std::complex<double>(0, 0) ||
				N(1, 2) != std::complex<double>(0, 0))
				N(2, 2) = 0, showeq = true;

			// Show equation
			if (showeq) {
				std::cout << std::endl;
				std::cout << "Nullify free temp if at least one linear coefficient is not a zero" << std::endl;
				outCurve(N);
				std::cout << std::endl;
			}

			// Reset flag
			showeq = false;

			// If all quadratic coefficients are zero and at least one linear coefficient is not a zero it is a line
			if (N(0, 0) == std::complex<double>(0, 0) &&
				N(1, 1) == std::complex<double>(0, 0) && (
				N(0, 2) != std::complex<double>(0, 0) ||
				N(1, 2) != std::complex<double>(0, 0))) {
				N(0, 2) = N(2, 0) = 0.5;
				N(1, 2) = N(2, 1) = 0;
				showeq = true;
			}

			// Show equation
			if (showeq) {
				std::cout << std::endl;
				std::cout << "If all quadratic coefficients are zero and" << std::endl;
				std::cout << "at least one linear coefficient is not a zero it is a line" << std::endl;
				outCurve(N);
				std::cout << std::endl;
			}

			// Reset flag
			showeq = false;

			// If all coefficients except for one quadratic are zero it is a line
			if (N(0, 2) == std::complex<double>(0, 0) &&
				N(1, 2) == std::complex<double>(0, 0) &&
				N(2, 2) == std::complex<double>(0, 0))
				for (uint32_t i = 0; i < 2; i++)
					if (N(i, i) != std::complex<double>(0, 0) &&
						N((i + 1) % 3, (i + 1) % 3) == std::complex<double>(0, 0)) {
						N(i, 2) = N(2, i) = 0.5;
						N(i, i) = 0;
						showeq = true;
					}

			// Show equation
			if (showeq) {
				std::cout << std::endl;
				std::cout << "If all coefficient except for one quadratic are zero it is a line" << std::endl;
				outCurve(N);
				std::cout << std::endl;
			}

			// Reset flag
			showeq = false;

			// Show equation
			std::cout << std::endl;
			std::cout << "Resulting equation" << std::endl;
			outCurve(N);
			std::cout << std::endl;

			// Quit current case
			break;

		// 2nd order surface exploration case
		case '2':

			// Introduce to user how to input 2nd order surface coefficients
			std::cout << "Enter the coefficients of the 2nd order surface as in the example:\n";
			std::cout << "If you have Ax^2 + By^2 + Cz^2 + Dxy + Exz + Fyz + Gx + Hy + Iz + J = 0\n";
			std::cout << "You should input <<A B C D E F G H I J>>\n";

			// Get 2nd order surface coefficients
			std::cin >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;

			// Divide repeating coefficients by 2
			D /= 2;
			E /= 2;
			F /= 2;
			G /= 2;
			H /= 2;
			I /= 2;

			// Construct coefficients matrix
			coeffs = MatrixXd(4, 4);
			coeffs <<
				A, D, E, G,
				D, B, F, H,
				E, F, C, I,
				G, H, I, J;

			// Calculate invariants
			subinvs = invariants(coeffs.block(0, 0, 3, 3));
			invs = invariants(coeffs);

			// Determine 2nd order surface by invariants combinations
			std::cout << std::endl;
			if (subinvs[2] != 0) {
				if (subinvs[1] > 0 && subinvs[0] * subinvs[2] > 0) {
					if (invs[3] < 0)
						std::cout << "Ellipsoid\n";
					else if (invs[3] > 0)
						std::cout << "Imaginary ellipsoid\n";
					else
						std::cout << "Point\n";
				} else if (subinvs[1] == 0 || subinvs[0] * subinvs[2] <= 0) {
					if (invs[3] > 0)
						std::cout << "One-sheet hyperboloid\n";
					else if (invs[3] < 0)
						std::cout << "Two-sheet hyperboloid\n";
					else
						std::cout << "Cone\n";
				}
			} else {
				if (invs[3] != 0) {
					if (invs[3] < 0)
						std::cout << "Elliptic paraboloid\n";
					else if (invs[3] > 0)
						std::cout << "Hyperbolic paraboloid\n";
				} else {
					if (subinvs[1] > 0) {
						if (subinvs[0] * invs[1] < 0)
							std::cout << "Elliptic cylinder\n";
						else if (subinvs[0] * invs[1] > 0)
							std::cout << "Imaginary elliptic cylinder\n";
						else if (invs[1] == 0)
							std::cout << "Line\n";
					} else if (subinvs[1] < 0) {
						if (invs[1] != 0)
							std::cout << "Hyperbolic cylinder\n";
						else
							std::cout << "Intersecting planes\n";
					} else {
						if (invs[1] != 0)
							std::cout << "Parabolic cylinder\n";
						else {
							if (invs[0] < 0)
								std::cout << "Parallel planes\n";
							else if (invs[0] > 0)
								std::cout << "Imaginary parallel planes\n";
							else
								std::cout << "Plane\n";
						}
					}
				}
			}

			// Calculate eigen vectors of quadratic coefficients submatrix
			es = EigenSolver<MatrixXd>(coeffs.block(0, 0, 3, 3));

			// Prepare eigen vectors matrix
			V = Matrix4cd::Identity();
			V.block(0, 0, 3, 3) = es.eigenvectors();

			// Apply new basis
			N = V.inverse() * coeffs * V;

			// Round near zero elements
			N = matRound(N);

			// Show equation
			std::cout << std::endl;
			std::cout << "Applying new basis" << std::endl;
			outSurface(N);
			std::cout << std::endl;

			// Completing the square
			for (uint32_t i = 0; i < 3; i++)
				if (N(i, i) != std::complex<double>(0, 0) &&
					N(i, 3) != std::complex<double>(0, 0)) {
					N(3, 3) -= N(i, 3) * N(i, 3) / N(i, i);
					N(i, 3) = N(3, i) = 0;
					showeq = true;
				}
			
			// Show equation
			if (showeq) {
				std::cout << std::endl;
				std::cout << "Completing squares" << std::endl;
				outSurface(N);
				std::cout << std::endl;
			}

			// Reset flag
			showeq = false;

			// Orthogonal substitution
			for (uint32_t i = 0; i < 3; i++)
				if (N((i + 1) % 3, 3) != std::complex<double>(0, 0) &&
					N((i + 2) % 3, 3) != std::complex<double>(0, 0)) {
					std::complex<double> u = sqrt(
						N((i + 1) % 3, 3) * N((i + 1) % 3, 3) +
						N((i + 2) % 3, 3) * N((i + 2) % 3, 3));
					N((i + 1) % 3, 3) = N(3, (i + 1) % 3) = u;
					N((i + 2) % 3, 3) = N(3, (i + 2) % 3) = 0;
					showeq = true;
				}

			// Show equation
			if (showeq) {
				std::cout << std::endl;
				std::cout << "Applying orthogonal substitution" << std::endl;
				outSurface(N);
				std::cout << std::endl;
			}

			// Reset flag
			showeq = false;

			// Nullify free temp if at least one linear coefficient is not a zero
			if (N(0, 3) != std::complex<double>(0, 0) ||
				N(1, 3) != std::complex<double>(0, 0) ||
				N(2, 3) != std::complex<double>(0, 0))
				N(3, 3) = 0, showeq = true;

			// Show equation
			if (showeq) {
				std::cout << std::endl;
				std::cout << "Nullify free temp if at least one linear coefficient is not a zero" << std::endl;
				outSurface(N);
				std::cout << std::endl;
			}

			// Reset flag
			showeq = false;

			// If all quadratic coefficients are zero and at least one linear coefficient is not a zero it is a plane
			if (N(0, 0) == std::complex<double>(0, 0) &&
				N(1, 1) == std::complex<double>(0, 0) &&
				N(2, 2) == std::complex<double>(0, 0) && (
				N(0, 3) != std::complex<double>(0, 0) ||
				N(1, 3) != std::complex<double>(0, 0) ||
				N(2, 3) != std::complex<double>(0, 0))) {
				N(0, 3) = N(3, 0) = 0.5;
				N(1, 3) = N(3, 1) = 0;
				N(2, 3) = N(3, 2) = 0;
				showeq = true;
			}

			// Show equation
			if (showeq) {
				std::cout << std::endl;
				std::cout << "If all quadratic coefficients are zero and" << std::endl;
				std::cout << "at least one linear coefficient is not a zero it is a plane" << std::endl;
				outSurface(N);
				std::cout << std::endl;
			}

			// Reset flag
			showeq = false;

			// If all coefficient except for one quadratic are zero it is a plane
			if (N(0, 3) == std::complex<double>(0, 0) &&
				N(1, 3) == std::complex<double>(0, 0) &&
				N(2, 3) == std::complex<double>(0, 0) &&
				N(3, 3) == std::complex<double>(0, 0))
				for (uint32_t i = 0; i < 3; i++)
					if (N(i, i) != std::complex<double>(0, 0) &&
						N((i + 1) % 3, (i + 1) % 3) == std::complex<double>(0, 0) &&
						N((i + 2) % 3, (i + 2) % 3) == std::complex<double>(0, 0)) {
						N(i, 3) = N(3, i) = 0.5;
						N(i, i) = 0;
						showeq = true;
					}

			// Show equation
			if (showeq) {
				std::cout << std::endl;
				std::cout << "If all coefficient except for one quadratic are zero it is a plane" << std::endl;
				outSurface(N);
				std::cout << std::endl;
			}

			// Reset flag
			showeq = false;

			// Show equation
			std::cout << std::endl;
			std::cout << "Resulting equation" << std::endl;
			outSurface(N);
			std::cout << std::endl;

			// Quit current case
			break;

		// Quit case
		case '0':

			// Be gentle
			std::cout << "\nGood bye\n";

			// Set running flag to false
			flag = false;

			// Quit current case
			break;

		// Default case
		default:

			// User's response is invalid
			std::cout << "\nYou chose the wrong answer. Try Again\n";

			// Quit current case
			break;
		}
	}
} /* End of 'CurvesAndSurfaces2ndOrder' function */

/* The main program function.
 * ARGUMENTS: None.
 * RETURNS: None.
 */
int main(void) {
	CurvesAndSurfaces2ndOrder();

	return 0;
} /* End of 'main' function */

/* END OF 'main.cpp' FILE */