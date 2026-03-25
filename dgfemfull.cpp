//DG-FEM Model: Suyog Vibhuti
//Uniform Advection Problem: dq/dt + c * dq/dx = 0, c is fluid velocity (constant)
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;
const int K = 32;
const int ORDER = 1; //so most things are 2 (basis functions)

void dgfem();
void forwardEuler(double (&a)[2][K], double aprime[2][K], double tStep);
void massMatrix(double (&mM)[2][2], bool inv);

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

	// Initial condition along sine wave, doing it this way is flawed since the lines between a points aren't really meant to be continuous, meant to be average of analytical solution
	/** double a[2][K];
    double tau = 2 * M_PI / static_cast<double>(K);
    // double length = 1; // determines length of system, how many periods basically
	for (int i = 0; i < K; i++) {
		a[0][i] = sin(tau * i);
		a[1][i] = sin(tau * (i + 1));
	} **/

    dgfem();
	
	// Exiting Program, Normal Operation Code
	return 0;
}

void dgfem() {
	double elementlist[K];
	double q[K];

	double aprime[2][K];
	double a[2][K];

	// Copied from main function for now!
	// Initial condition along sine wave, doing it this way is flawed since the lines between a points aren't really meant to be continuous, meant to be average of analytical solution
    double tau = 2 * M_PI / static_cast<double>(K);
	for (int i = 0; i < K; i++) {
		a[0][i] = sin(tau * static_cast<double>(i));
		a[1][i] = sin(tau * static_cast<double>(i + 1));
	}

    // First part of output, for original values
    ofstream file("results.txt");
    for (int i = 0; i < K; i++) {
        file << a[0][i] << "," << a[1][i] << "\n";
    }
    file << "\n";

	// aprime = invM[cKa - f]

	int c = 3; // dummy value for fluid velocity, make dgfem parameter later
	for (int i = 1; i < K; i++) {
		// Formulas derived from page in notebook, lots of matrix multiplication
		aprime[0][i] = c * (-3 * a[0][i] - a[1][i] + 4 * a[1][i - 1]);
		aprime[1][i] = c * (3 * a[0][i] - a[1][i] - 2 * a[1][i - 1]);
	}
	// Wraparound condition, so that values on right affect values on left
	aprime[0][0] = c * (-3 * a[0][0] - a[1][0] + 4 * a[1][K - 1]);
	aprime[1][0] = c * (3 * a[0][0] - a[1][0] - 2 * a[1][K - 1]);

	// To-Do: Develop ODE integrator (forward euler or RK4), apply to aprime values until desired time t, make q list: q_i = a_0i(1 - x) + a_1i(x). Output q list?
    double time = 1;
    double tStep = 0.01;
    for (int secondsCount = 0; secondsCount < time; secondsCount++) {
        for (int count = 0; count < (1.0 / tStep); count++) {
			cout << a[0][6] << "\n";
            forwardEuler(a, aprime, tStep);
            for (int i = 1; i < K; i++) {
		        // Formulas derived from page in notebook, lots of matrix multiplication
		        aprime[0][i] = c * (-3 * a[0][i] - a[1][i] + 4 * a[1][i - 1]);
		        aprime[1][i] = c * (3 * a[0][i] - a[1][i] - 2 * a[1][i - 1]);
	        }
	        // Wraparound condition, so that values on right affect values on left
	        aprime[0][0] = c * (-3 * a[0][0] - a[1][0] + 4 * a[1][K - 1]);
	        aprime[1][0] = c * (3 * a[0][0] - a[1][0] - 2 * a[1][K - 1]);

            // write results into file
            for (int i = 0; i < K; i++) {
                file << a[0][i] << "," << a[1][i] << "\n";
            }
            file << "\n";
        }

        /** // write results into file
        for (int i = 0; i < K; i++) {
            file << a[0][i] << "," << a[1][i] << "\n";
        }
        file << "\n"; **/
    }

    file.close();
}

void forwardEuler(double (&a)[2][K], double aprime[2][K], double tStep) {
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < K; j++) {
            a[i][j] = a[i][j] + tStep * aprime[i][j];
        }
    }
}

void massMatrix(double (&mM)[2][2], bool inv) {
	if (inv) {
		// For returning the inverse mass matrix
		mM[0][0] = 4;
		mM[1][1] = 4;
		mM[0][1] = -2;
		mM[1][0] = -2;
	} else {
		// For returning the standard mass matrix
		mM[0][0] = 1/3;
		mM[1][1] = 1/3;
		mM[0][1] = 1/6;
		mM[1][0] = 1/6;
	}
}

void stiffnessMatrix(double (&sM)[2][2]) {
	// For simplified form, order 1
	// first term is row, second term is column
	sM[0][0] = -1/2;
	sM[0][1] = -1/2;
	sM[1][0] = 1/2;
	sM[1][1] = 1/2;
}

int fluxTerm() {
	double fT[2];
	// For simplified form, order 1
	// fT[0] = c * a1;
	// fT[1] = c * a2;
	return 0;
}