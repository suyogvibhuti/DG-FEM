//DG-FEM Model: Suyog Vibhuti
//Uniform Advection Problem: dq/dt + c * dq/dx = 0, c is fluid velocity (constant)
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;
const int K = 64;
const int ORDER = 2; // Moving to arbitrary order, generalizing matrices and procedures

void dgfem(double fluidVelocity, double length);

int main() {
    // For arbitrary order need to use lagrange interpolation on nodal system (at node pts approximation equals analytical solution)
    // Runge's phenomenon for equally spaced nodes, use legendre interpolation points instead for more points at edges
    // Legendre interpolation: take a set of legendre polynomials of given order, roots of polynomials form set of interpolation pts
    // For HO matrices (mass, stiffness), in order to integrate arbitrary functions use quadrature (creating a new interpolation of function using lagrange interpolation)
    // Try hermite interpolation quadrature (2M - 1 order), doesnt need value for first deriv. w/ legendre-lagrange bases
    // hermite interpolation: integral of f ~ sum from 0 to M of f(xi)*wi

    // Running dgfem function, dummy value for fluid velocity
    dgfem(2.0, 32.0);

    // Exiting program, normal operation code
    return 0;
}

void dgfem(double fluidVelocity, double length) {
    // designed to solve for reference element, then to apply to elements in region using jacobians
    const int Nv = K + 1; // number of vertices
    double vx[Nv]; // row vector, each value represents endpoint of element on region
    int EToV[K][2]; // each row represents element: contains two values from vx for left and right points
    
    // affine mapping (for element k): x(r) = xl + ((1 + r)/2) * h, h = (xr - xl)
        // r is between -1 and 1
    
    int FToV[2 * K][Nv];
    int FToV_T[Nv][2 * K]; // Transposed FToV
    for (int i = 0; i < K; i++) {
        FToV[2 * i][EToV[i][0] - 1], FToV_T[EToV[i][0] - 1][2 * i] = 1;
        FToV[2 * i + 1][EToV[i][1] - 1], FToV_T[EToV[i][1] - 1][2 * i + 1] = 1;
    }

    // multiplying FToV by FToV_T and subtracting the identity matrix
    int FToF[2 * K][2 * K];
    for (int i = 0; i < 2 * K; i++) {
        for (int j = 0; j < 2 * K; j++) {
            for (int k = 0; k < Nv; k++) {
                FToF[i][j] += FToV[i][k] * FToV_T[k][j];
            }
        }
        FToF[i][i] = 0; // effectively works as identity matrix subtraction operation 
    }

    cout << "hello!";
}

int* JacobiGQ(int alpha, int beta, int N) {
    if (N == 0) {
        return {0};
    }
    int* h1 = new int[N + 1];
    for (int n = 0; n < N + 1; n++) {
        h1[n] = 2 * n + alpha + beta;
    }

    double** J = new double*[N + 1];
    for (int i = 0; i < N + 1; i++) {
        J[i] = new double[N + 1];
        J[i][i] = -static_cast<double>(alpha * alpha - beta * beta) / static_cast<double>((h1[i] + 2) * h1[i]); // For a=b=0, equals 0. for gauss-lobatto, since a=b=1, doesn't equal 0.
    }
    for (int i = 0; i < N; i++) {
        J[i][i + 1], J[i + 1][i] = (2.0 / static_cast<double>(h1[i])) * sqrt(static_cast<double>(i * (i + alpha + beta) * (i + alpha) * (i + beta)) / static_cast<double>((h1[i] - 1) * (h1[i] + 1)));
    }
    if (alpha + beta == 0) {
        J[0][0] = 0;
    }
}