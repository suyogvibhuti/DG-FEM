//DG-FEM Model: Suyog Vibhuti
//Uniform Advection Problem: dq/dt + c * dq/dx = 0, c is fluid velocity (constant)
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;
const int K = 32;
const int ORDER = 1; // so most things are 2 (basis functions)

void dgfem(double fluidVelocity);
void forwardEuler(double (&a)[2][K], double aprime[2][K], double tStep);
void massMatrix(double (&mM)[2][2], bool inv);
void stiffnessMatrix(double (&sM)[2][2]);
void fluxTerm(double (&fT)[2], double c, double a1, double a2);
void matrixColMult(double (&matrix2D)[2][2], double (&matrix1D)[2], double (&resultMatrix)[2]);
void numericalIntegration(double a[2][K], double c, int i, double (&resultMatrix)[2]);

int main() {
	// differential cell dx by dy, possesses dq/dt, flux on each side fx(q) for left and fx'(q) for right
	// dq/dt handled with ODE solver, df(q)/dx handled with PDE solver
	// flux is linear, f(q) = cq, df(q)/dx = c* dq/dx
	// integral domain over element xl to xr
	// flux at boundary ^f(q+, q-) = cq- (q- is left side), cq+ if c is negative
	// approx. solution represented as ~q(x) = sum i to M of a_i * psi_i(x)
		// M = 1 for simplified basis vectors of order 1, ramp functions
		// psi_0(x) = 1 - x, psi_1(x) = x, ~q(x) = a0(1-x) + a1(x), line between a0 and a1
			// ql = a0, qr = a1, benefit of discontinuity is getting closer to analytical solution
	// inner product of analytical solution with basis should equal inner product of approx. sol. with basis
	// instead of dealing with df/dx, use integration by parts to deal with ds/dx, s is test function
	

	// Things needed for each element k within set of elements K
		// position: x in k, xl <= x <= xr

	// Boundary condition: q(0, t) = g(t), where g(t) is arbitrary function

	// Running dgfem function, dummy value of 3 for fluid velocity
    dgfem(3);
	
	// Exiting Program, Normal Operation Code
	return 0;
}

void dgfem(double fluidVelocity) {
	double elementlist[K];
	double q[K];

	double aprime[2][K];
	double a[2][K];

	// Initial condition along sine wave, doing it this way is flawed since the lines between a points aren't really meant to be continuous, meant to be average of analytical solution
    double tau = 2 * M_PI / static_cast<double>(K);
	for (int i = 0; i < K; i++) {
		a[0][i] = sin(tau * static_cast<double>(i));
		a[1][i] = sin(tau * static_cast<double>(i + 1));
		// Changed initial condition to square wave, doubled frequency
		/** a[0][i] = -1;
		a[1][i] = -1;
		if (sin(2 * tau * static_cast<double>(i)) > 0) {
			a[0][i] = 1;
		}
		if (sin(2 * tau * static_cast<double>(i + 1)) > 0) {
			a[1][i] = 1;
		} **/
	}

    // First part of output, for original values
    ofstream file("results.txt");
	file << K << "\n" << "\n"; // For graphing script, communicated number of elements to graph
    for (int i = 0; i < K; i++) {
        file << a[0][i] << "," << a[1][i] << "\n";
    }
    file << "\n";

	double c = fluidVelocity;
	double aprimei[2];
	for (int i = 0; i < K; i++) {
		// Using numerical integrator function for formula aprime = invM[cKa - f]. Includes wraparound condition.
		numericalIntegration(a, c, i, aprimei);
		aprime[0][i] = aprimei[0];
		aprime[1][i] = aprimei[1];
	}

	// To-Do: Develop ODE integrator (forward euler or RK4), apply to aprime values until desired time t
    double time = 100;
    double tStep = 0.01;
	int writeStep = 1;
	double deltaX = 6.28 / K;
	double CFL = c * tStep / deltaX;
	int elapsedTimeCounter = 1;
	int writeCount = 1;
    for (int secondsCount = 0; secondsCount < time; secondsCount++) {
        for (int count = 0; count < (1.0 / tStep); count++) {
			// cout << a[0][6] << "\n";
            forwardEuler(a, aprime, tStep);
            for (int i = 0; i < K; i++) {
		        // Using numerical integrator function for formula aprime = invM[cKa - f]. Includes wraparound condition.
				numericalIntegration(a, c, i, aprimei);
				aprime[0][i] = aprimei[0];
				aprime[1][i] = aprimei[1];
	        }

			// write results into file
			if (elapsedTimeCounter % (writeStep * static_cast<int>(1.0 / tStep)) == 0) {
				for (int i = 0; i < K; i++) {
                	file << a[0][i] << "," << a[1][i] << "\n";
            	}
            	file << "\n";
				writeCount++;
				cout << a[0][6] << "\n";
			}
			elapsedTimeCounter++;
        }

		/** // write results into file
		if (secondsCount % writeStep == 0) {
			for (int i = 0; i < K; i++) {
                file << a[0][i] << "," << a[1][i] << "\n";
            }
            file << "\n";
		} **/
    }

    file.close();
	cout << writeCount;
}

void forwardEuler(double (&a)[2][K], double aprime[2][K], double tStep) {
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < K; j++) {
            a[i][j] = a[i][j] + tStep * aprime[i][j];
        }
    }
}

void sspRK3(double (&a)[2][K], double aprime[2][K], double tStep) {
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < K; j++) {
			double f1 = 0.5 * a[i][j] + 0.5 * (a[i][j] + tStep * aprime[i][j]);
			double f2 = 0.5 * f1 + 0.5 * (f1 + tStep * aprime[i][j]);
		}
	}
}

void massMatrix(double (&mM)[2][2], bool inv) {
	if (inv) {
		// For returning the inverse mass matrix
		mM[0][0] = 4.0;
		mM[1][1] = 4.0;
		mM[0][1] = -2.0;
		mM[1][0] = -2.0;
	} else {
		// For returning the standard mass matrix
		mM[0][0] = 1.0 / 3.0;
		mM[1][1] = 1.0 / 3.0;
		mM[0][1] = 1.0 / 6.0;
		mM[1][0] = 1.0 /6.0;
	}
}

void stiffnessMatrix(double (&sM)[2][2]) {
	// For simplified form, order 1
	// first term is row, second term is column
	sM[0][0] = -1.0 / 2.0;
	sM[0][1] = -1.0 / 2.0;
	sM[1][0] = 1.0 / 2.0;
	sM[1][1] = 1.0 / 2.0;
}

void fluxTerm(double (&fT)[2], double c, double a1, double a2) {
	// For simplified form, order 1
	fT[0] = -1 * c * a1;
	fT[1] = c * a2;
}

void matrixColMult(double (&matrix2D)[2][2], double (&matrix1D)[2], double (&resultMatrix)[2]) {
	// For simplified form, order 1. To extend, don't hardcode the size of the matrices and use loops
	resultMatrix[0] = matrix1D[0] * matrix2D[0][0] + matrix1D[1] * matrix2D[0][1];
	resultMatrix[1] = matrix1D[0] * matrix2D[1][0] + matrix1D[1] * matrix2D[1][1];
}

void numericalIntegration(double a[2][K], double c, int i, double (&resultMatrix)[2]) {
	// aprime = invM[cKa - f], K refers to stiffness matrix sM
	double sM[2][2];
	stiffnessMatrix(sM);
	double acol[2] = {a[0][i], a[1][i]};
	double Ka[2];
	matrixColMult(sM, acol, Ka);
	double f[2];
	// For wraparound condition, so that values on right affect values on left
	if (i > 0) {
		fluxTerm(f, c, a[1][i - 1], a[1][i]);
	} else {
		fluxTerm(f, c, a[1][K - 1], a[1][i]);
	}
	double term2[2] = {c * Ka[0] - f[0], c * Ka[1] - f[1]};
	double invM[2][2];
	massMatrix(invM, true);
	matrixColMult(invM, term2, resultMatrix);
}