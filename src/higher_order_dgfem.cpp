// DG-FEM Model: Suyog Vibhuti
// Uniform Advection Problem: dq/dt + c * dq/dx = 0, c is fluid velocity (constant)
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace std;
using namespace Eigen;
const int K = 64;
const int ORDER = 2; // Moving to arbitrary order, generalizing matrices and procedures

void dgfem(double fluidVelocity, double length);
void startup(VectorXd VX, MatrixXd EToV, int K);
VectorXd JacobiGQ(int alpha, int beta, int N);
VectorXd JacobiGL(int alpha, int beta, int N);
VectorXd JacobiP(VectorXd x, int x_length, int alpha, int beta, int N);
MatrixXd Vandermonde1D(int N, VectorXd r, int r_length);
VectorXd GradJacobiP(VectorXd r, int r_length, int alpha, int beta, int N);
MatrixXd GradVandermonde1D(VectorXd r, int r_length, int N);
MatrixXd Dmatrix1D(VectorXd r, int r_length, int N, MatrixXd Vandermonde);
MatrixXd Lift1D(int r_length, int N, int Nfaces, int Nfp, MatrixXd Vandermonde);
MatrixXd Jacobian1D(VectorXd x, MatrixXd Dr, int r_length);

int main() {
    // For arbitrary order need to use lagrange interpolation on nodal system (at node pts approximation equals analytical solution)
    // Runge's phenomenon for equally spaced nodes, use legendre interpolation points instead for more points at edges
    // Legendre interpolation: take a set of legendre polynomials of given order, roots of polynomials form set of interpolation pts
    // For HO matrices (mass, stiffness), in order to integrate arbitrary functions use quadrature (creating a new interpolation of function using lagrange interpolation)
    // Try hermite interpolation quadrature (2M - 1 order), doesnt need value for first deriv. w/ legendre-lagrange bases
    // hermite interpolation: integral of f ~ sum from 0 to M of f(xi)*wi

    // Running dgfem function, dummy value for fluid velocity
    int K = 3;
    VectorXd VX(K + 1);
    MatrixXd EToV(K, 2);
    startup(VX, EToV, K);
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

    cout << "hello!" << "\n";
}

void startup(VectorXd VX, MatrixXd EToV, int K) {
    // Setup script to define operators, grid, and connection faces

    // Constants
    int N = 16;
    int numFaces = 2;
    int Nfp = 1;

    // Basic LGL grid for reference element r
    VectorXd r(N + 1);
    int r_length = N + 1;
    r = JacobiGL(0, 0, N);
    cout << "<";
    for (int i = 0; i < N; i++) {
        cout << r(i) << ", ";
    }
    cout << r(N) << ">\n";

    // Reference element matrices
    MatrixXd V(N + 1, N);
    V = Vandermonde1D(N, r, r_length);
    MatrixXd invV(N, N + 1);
    invV = V.inverse();

    MatrixXd Dr(r_length, r_length);
    Dr = Dmatrix1D(r, r_length, N, V);

    // Surface integral terms
    MatrixXd LIFT(r_length, numFaces * Nfp);
    LIFT = Lift1D(r_length, N, numFaces, Nfp, V);

    VectorXd va(K);
    VectorXd vb(K);
    va = EToV.col(0).array() - 1;
    vb = EToV.col(1).array() - 1;
    VectorXd x(N + 1);
    x = VectorXd::Ones(N + 1) * VX(Eigen::placeholders::all, va);
    VectorXd r_intermediary(r_length);
    r_intermediary = r.array() + 1;
    r_intermediary = 0.5 * r_intermediary * (VX(Eigen::placeholders::all, vb) - VX(Eigen::placeholders::all, va));
    x = x.array() + r_intermediary.array();

    VectorXd J(r_length);
    J = Jacobian1D(x, Dr, r_length);
    VectorXd rx(r_length);
    for (int i = 0; i < r_length; i++) {
        rx(i) = 1.0 / J(i);
    }
}

VectorXd JacobiGQ(int alpha, int beta, int N) {
    if (N == 0) {
        VectorXd x(1);
        x(0) = static_cast<double>(alpha - beta) / static_cast<double>(alpha + beta + 2);
        return x;
    }
    int* h1 = new int[N + 1];
    for (int n = 0; n < N + 1; n++) {
        h1[n] = 2 * n + alpha + beta;
    }

    double** J = new double*[N + 1];
    for (int i = 0; i < N + 1; i++) {
        J[i] = new double[N + 1];
        J[i][i] = -static_cast<double>(alpha * alpha - beta * beta) / ((h1[i] + 2) * h1[i]); // For a=b=0, equals 0. for gauss-lobatto, since a=b=1, doesn't equal 0.
    }
    for (int i = 1; i < N + 1; i++) {
        J[i - 1][i] = (2.0 / static_cast<double>(h1[i])) * sqrt(static_cast<double>(i * (i + alpha + beta) * (i + alpha) * (i + beta)) / static_cast<double>((h1[i] - 1) * (h1[i] + 1)));
        J[i][i - 1] = (2.0 / static_cast<double>(h1[i])) * sqrt(static_cast<double>(i * (i + alpha + beta) * (i + alpha) * (i + beta)) / static_cast<double>((h1[i] - 1) * (h1[i] + 1)));
    }
    if (alpha + beta == 0) {
        J[0][0] = 0;
    }

    MatrixXd mat(N + 1, N + 1);
    for (int i = 0; i < N + 1; i++) {
        for (int j = 0; j < N + 1; j++) {
            // cout << J[i][j] << " ";
            mat(i, j) = J[i][j];
        }
        // cout << "\n";
    }
    // cout << "\n";
    // cout << mat;
    EigenSolver<MatrixXd> solver(mat);
    VectorXd x(N + 1);
    VectorXcd J_eigenvalues = solver.eigenvalues();
    for (int i = 0; i < N + 1; i++) {
        x(i) = J_eigenvalues(i).real();
    }
    cout << solver.eigenvalues().real() << "\n";
    return x;
}

VectorXd JacobiGL(int alpha, int beta, int N) {
    VectorXd x(N + 1);
    x(0) = -1.0;
    x(N) = 1.0;
    if (N > 1) {
        VectorXd xint(N - 1);
        xint = JacobiGQ(alpha + 1, beta + 1, N - 2);
        for (int i = 1; i < N; i++) {
            x(i) = xint(i - 1);
        }
    }
    return x;
}

VectorXd JacobiP(VectorXd x, int x_length, int alpha, int beta, int N) {
    MatrixXd PL(N + 1, x_length);

    double gamma0 = pow(2.0, alpha + beta + 1) * tgamma(alpha + 1) * tgamma(beta + 1) / tgamma(alpha + beta + 2);
    double PL_0 = 1.0 / sqrt(gamma0);
    for (int i = 0; i < x_length; i++) {
        PL(0, i) = PL_0;
    }
    if (N == 0) {
        VectorXd PL_row0(x_length);
        for (int i = 0; i < x_length; i++) {
            PL_row0(i) = PL(0, i);
        }
        return PL_row0;
    }

    double gamma1 = gamma0 * static_cast<double>((alpha + 1) * (beta + 1)) / static_cast<double>((alpha + beta + 3));
    for (int i = 0; i < x_length; i++) {
        PL(1, i) = (static_cast<double>(alpha + beta + 2) * x(i) + static_cast<double>(alpha- beta)) / (2.0 * sqrt(gamma1));
    }
    if (N == 1) {
        VectorXd PL_row1(x_length);
        for (int i = 0; i < x_length; i++) {
            PL_row1(i) = PL(1, i);
        }
        return PL_row1;
    }

    double aold = 2.0 * sqrt(static_cast<double>((alpha + 1) * (beta + 1)) / static_cast<double>(alpha + beta + 3)) / (2 + alpha + beta); // aold = a1
    for (int i = 0; i < N - 1; i++) {
        int h1 = 2 * i + alpha + beta;
        double anew = 2.0 * sqrt(static_cast<double>((i + 1) * (i + 1 + alpha + beta) * (i + 1 + alpha) * (i + 1 + beta)) / static_cast<double>((h1 + 1) * (h1 + 3))) / (h1 + 2);
        double bnew = -static_cast<double>(alpha * alpha - beta * beta) / (h1 * (h1 + 2));
        for (int j = 0; j < x_length; j++) {
            PL(i + 2, j) = (-aold * PL(i, j) + (x(j) - bnew) * PL(i + 1, j)) / anew;
        }
        aold = anew;
    }

    VectorXd PL_rowN(x_length);
    for (int i = 0; i < x_length; i++) {
        PL_rowN(i) = PL(N, i);
    }
    return PL_rowN;
}

MatrixXd Vandermonde1D(int N, VectorXd r, int r_length) {
    // Sets up Vandermonde matrix, used to convert between nodal and modal formulations
    MatrixXd V1D(r_length, N + 1);

    VectorXd V1Dcol(r_length);
    for (int j = 0; j < N + 1; j++) {
        V1Dcol = JacobiP(r, r_length, 0, 0, j);
        for (int i = 0; i < r_length; i++) {
            V1D(i, j) = V1Dcol(i);
        }
    }

    return V1D;
}

VectorXd GradJacobiP(VectorXd r, int r_length, int alpha, int beta, int N) {
    // Finds derivative of jacobi polynomial based on identity
    VectorXd deltaP(r_length);
    if (!(N == 0)) {
        deltaP = JacobiP(r, r_length, alpha + 1, beta + 1, N - 1);
        for (int i = 0; i < r_length; i++) {
            deltaP(i) = deltaP(i) * sqrt(N * (N + alpha + beta + 1));
        }
    }

    return deltaP;
}

MatrixXd GradVandermonde1D(VectorXd r, int r_length, int N) {
    // Compiles Jacobi derivatives to create Vandermonde gradient
    MatrixXd dVr(r_length, N + 1);

    for (int i = 0; i < N + 1; i++) {
        VectorXd dVrCol(r_length);
        dVrCol = GradJacobiP(r, r_length, 0, 0, i);
        for (int j = 0; j < r_length; j++) {
            dVr(j, i) = dVrCol(j);
        }
    }

    return dVr;
}

MatrixXd Dmatrix1D(VectorXd r, int r_length, int N, MatrixXd Vandermonde) {
    // Finds differentiation matrix for which M * D = S
    MatrixXd Vr(r_length, N + 1);
    Vr = GradVandermonde1D(r, r_length, N);

    MatrixXd Dmatrix(r_length, r_length);
    Dmatrix = Vr * Vandermonde.inverse();

    return Dmatrix;
}

MatrixXd Lift1D(int r_length, int N, int Nfaces, int Nfp, MatrixXd Vandermonde) {
    MatrixXd EMat(N + 1, Nfaces * Nfp);
    EMat(0, 0) = 1.0;
    EMat(N, 1) = 1.0;

    MatrixXd result_matrix(N + 1, Nfaces * Nfp);
    result_matrix = Vandermonde.transpose() * EMat;
    MatrixXd result_matrix_2(r_length, Nfaces * Nfp);
    result_matrix_2 = Vandermonde * result_matrix;

    return result_matrix_2;
}

MatrixXd Jacobian1D(VectorXd x, MatrixXd Dr, int r_length) {
   // Used for affine mapping from x to r, rx calculation put in startup
   
   // Calculating Jacobian
   VectorXd J(r_length);
   J = Dr * x;

   return J;
}