//DG-FEM Model: Suyog Vibhuti
//Uniform Advection Problem: dq/dt + c * dq/dx = 0, c is fluid velocity (constant)
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;
const int K = 64;
const int ORDER = 1; // so most things are 2 (basis functions)

void dgfem(double fluidVelocity, double length);
void forwardEuler(double (&a)[2][K + 1], double aprime[2][K], double tStep);
void sspRK3(double (&a)[2][K + 1], double aprime[2][K], double tStep);
void massMatrix(double (&mM)[2][2], bool inv);
void stiffnessMatrix(double (&sM)[2][2]);
void fluxTerm(double (&fT)[2], double c, double a1, double a2);
void matrixColMult(double (&matrix2D)[2][2], double (&matrix1D)[2], double (&resultMatrix)[2]);
void numericalIntegration(double a[2][K + 1], double c, int i, double (&resultMatrix)[2]);

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

	// Running dgfem function, dummy value for fluid velocity
    dgfem(2.0, 32.0);
	
	// Exiting Program, Normal Operation Code
	return 0;
}

void dgfem(double fluidVelocity, double length) {
	double elementlist[K];
	double q[K];

	double aprime[2][K];
	double a[2][K + 1]; // Trying to create dummy values on either side so that direction of upwind flux is irrelevant

	// Initial condition along sine wave, doing it this way is flawed since the lines between a points aren't really meant to be continuous, meant to be average of analytical solution
    double tau = 2 * M_PI / static_cast<double>(K);
	for (int i = 0; i < K; i++) {
		/** a[0][i] = sin(tau * static_cast<double>(i));
		a[1][i] = sin(tau * static_cast<double>(i + 1)); **/
		// Changed initial condition to square wave, doubled frequency
		a[1][i] = -1;
		a[0][i + 1] = -1;
		if (sin(2 * tau * static_cast<double>(i)) > 0) {
			a[1][i] = 1;
		}
		if (sin(2 * tau * static_cast<double>(i + 1)) > 0) {
			a[0][i + 1] = 1;
		}
	}

	// Creating copy values on either side
	a[0][0] = a[0][K];
	a[1][K] = a[1][0];

    // First part of output, for original values
    ofstream file("results.txt");
	file << K << "\n" << "\n"; // For graphing script, communicated number of elements to graph
    for (int i = 0; i < K; i++) {
        file << a[1][i] << "," << a[0][i + 1] << "\n";
    }
    file << "\n";

	double c = fluidVelocity * static_cast<double>(K) / length; // Every element has a length of "1" if taken literally, this enforces true length
	// fluidVelocity is how fast the fluid is traveling in m/s, c is measured in elements/s, length/K is the length of an element
	cout << c << " is element velocity\n";
	double aprimei[2];
	for (int i = 0; i < K; i++) {
		// Using numerical integrator function for formula aprime = invM[cKa - f]. Includes wraparound condition.
		numericalIntegration(a, c, i, aprimei);
		aprime[0][i] = aprimei[0];
		aprime[1][i] = aprimei[1];
	}

	// To-Do: Develop ODE integrator (forward euler or RK4), apply to aprime values until desired time t
    double time = 100;
    double tStep = 0.00005;
	int writeStep = 1;
	double deltaX = length / K;
	double CFL = fluidVelocity * tStep / deltaX; // CFL = c * tStep as well, CFL <= 1 for model to work properly
	cout << CFL << " is CFL\n";
	int elapsedTimeCounter = 1;
	int writeCount = 1;
    for (int secondsCount = 0; secondsCount < time; secondsCount++) {
        for (int count = 0; count < (1.0 / tStep); count++) {
			// cout << a[0][6] << "\n";
            forwardEuler(a, aprime, tStep);
			// sspRK3(a, aprime, tStep);

			// reset dummy values on ends of a
			a[0][0] = a[0][K];
			a[1][K] = a[1][0];

            for (int i = 0; i < K; i++) {
		        // Using numerical integrator function for formula aprime = invM[cKa - f]. Includes wraparound condition.
				numericalIntegration(a, c, i, aprimei);
				aprime[0][i] = aprimei[0];
				aprime[1][i] = aprimei[1];
	        }

			// write results into file
			if (elapsedTimeCounter % (writeStep * static_cast<int>(1.0 / tStep)) == 0) {
				for (int i = 0; i < K; i++) {
                	file << a[1][i] << "," << a[0][i + 1] << "\n";
            	}
            	file << "\n";
				writeCount++;
				cout << a[1][6] << "\n";
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

void forwardEuler(double (&a)[2][K + 1], double aprime[2][K], double tStep) {
    for (int i = 0; i < 2; i++) {
		int shifted_i = 1 - i; // in order to ensure aprime properly maps onto a
        for (int j = 0; j < K; j++) {
            a[shifted_i][j + i] = a[shifted_i][j + i] + tStep * aprime[i][j];
        }
    }
}

void sspRK3(double (&a)[2][K + 1], double aprime[2][K], double tStep) {
	// Basically identical to forward euler for order of 1, since it just repeats using aprime
	// Probably needs numerical integration updates each step of the way
	for (int i = 0; i < 2; i++) {
		int shifted_i = 1 - i; // in order to ensure aprime properly maps onto a
		for (int j = 0; j < K; j++) {
			double f1 = a[shifted_i][j + i] + tStep * aprime[i][j];
			double f2 = 0.75 * a[shifted_i][j + i] + 0.25 * (f1 + (tStep * aprime[i][j]));
			a[shifted_i][j + i] = (1.0 / 3.0) * a[shifted_i][j + i] + (2.0 / 3.0) * (f2 + (tStep * aprime[i][j]));
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

void numericalIntegration(double a[2][K + 1], double c, int i, double (&resultMatrix)[2]) {
	// aprime = invM[cKa - f], K refers to stiffness matrix sM
	double sM[2][2];
	stiffnessMatrix(sM);
	double acol[2] = {a[1][i], a[0][i + 1]};
	double Ka[2];
	matrixColMult(sM, acol, Ka);
	double f[2];
	// No need for wraparound condition anymore
	fluxTerm(f, c, a[0][i], a[0][i + 1]);
	double term2[2] = {c * Ka[0] - f[0], c * Ka[1] - f[1]};
	double invM[2][2];
	massMatrix(invM, true);
	matrixColMult(invM, term2, resultMatrix);
}