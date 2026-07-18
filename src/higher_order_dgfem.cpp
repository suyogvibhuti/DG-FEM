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

// Constants
const int K = 5;
const int N = 8;
const int numFaces = 2;
const int Nfp = 1;
const double node_tolerance = pow(10.0, -10.0);
const int ORDER = 2; // Moving to arbitrary order, generalizing matrices and procedures

// Global Variables
VectorXd r(N + 1);
int r_length = N + 1;
MatrixXd V(N + 1, N);
MatrixXd invV(N, N + 1);
MatrixXd Dr(r_length, r_length);
MatrixXd LIFT(r_length, numFaces * Nfp);
MatrixXd x(N + 1, K);
MatrixXd J(r_length, K);
MatrixXd rx(r_length, K);

// Low-storage Runge-Kutta coeffs
VectorXd rk4a(5);
VectorXd rk4b(5);
VectorXd rk4c(5);

// void dgfem(double fluidVelocity, double length);
void startup(VectorXd VX, MatrixXd EToV, int K);
VectorXd JacobiGQ(int alpha, int beta, int N);
VectorXd JacobiGL(int alpha, int beta, int N);
VectorXd JacobiP(VectorXd x, int x_length, int alpha, int beta, int N);
MatrixXd Vandermonde1D(int N, VectorXd r, int r_length);
VectorXd GradJacobiP(VectorXd r, int r_length, int alpha, int beta, int N);
MatrixXd GradVandermonde1D(VectorXd r, int r_length, int N);
MatrixXd Dmatrix1D(VectorXd r, int r_length, int N, MatrixXd Vandermonde);
MatrixXd Lift1D(int r_length, int N, int Nfaces, int Nfp, MatrixXd Vandermonde);
MatrixXd Jacobian1D(MatrixXd x, MatrixXd Dr, int r_length);
MatrixXd Normals1D(int Nfp, int Nfaces, int K);
void AdvecRHS1D(VectorXd u, double time, double a);
MatrixXd Advec1D(MatrixXd u, double finalTime);

int main() {
    // For arbitrary order need to use lagrange interpolation on nodal system (at node pts approximation equals analytical solution)
    // Runge's phenomenon for equally spaced nodes, use legendre interpolation points instead for more points at edges
    // Legendre interpolation: take a set of legendre polynomials of given order, roots of polynomials form set of interpolation pts
    // For HO matrices (mass, stiffness), in order to integrate arbitrary functions use quadrature (creating a new interpolation of function using lagrange interpolation)
    // Try hermite interpolation quadrature (2M - 1 order), doesnt need value for first deriv. w/ legendre-lagrange bases
    // hermite interpolation: integral of f ~ sum from 0 to M of f(xi)*wi

    VectorXd VX(K + 1);
    MatrixXd EToV(K, 2);

    // Uniform grid generation case
    double xMax = 32.0;
    double xMin = 0.0;
    for (int i = 0; i < K + 1; i++) {
        VX(i) = xMin + (xMax - xMin) * static_cast<double>(i) / static_cast<double>(K);
    }
    for (int i = 0; i < K; i++) {
        EToV(i, 0) = i;
        EToV(i, 1) = i + 1;
    }

    startup(VX, EToV, K);

    // Setting initial conditions for fluid velocity
    MatrixXd u(N + 1, K);
    for (int i = 0; i < N + 1; i++) {
        for (int j = 0; j < K; j++) {
            u(i, j) = sin(x(i, j));
        }
    }

    // Setting Runge-Kutta coeffs
    rk4a << 0.0, -567301805773.0/1357537059087.0, -2404267990393.0/2016746695238.0, -3550918686646.0/2091501179385.0, -1275806237668.0/842570457699.0;
    rk4b << 1432997174477.0/9575080441755.0, 5161836677717.0/13612068292357.0, 1720146321549.0/2090206949498.0, 3134564353537.0/4481467310338.0, 2277821191437.0/14882151754819.0;
    rk4c << 0.0, 1432997174477.0/9575080441755.0, 2526269341429.0/6820363962896.0, 2006345519317.0/3224310063776.0, 2802321613138.0/2924317926251.0;

    // Solve problem and print result
    double finalTime = 10;
    u = Advec1D(u, finalTime);
    cout << u;

    // Exiting program, normal operation code
    return 0;
}

/* void dgfem(double fluidVelocity, double length) {
    // designed to solve for reference element, then to apply to elements in region using jacobians
    const int Nv = K + 1; // number of vertices
    double vx[Nv]; // row vector, each value represents endpoint of element on region
    int EToV[K][2]; // each row represents element: contains two values from vx for left and right points
    
    // affine mapping (for element k): x(r) = xl + ((1 + r)/2) * h, h = (xr - xl)
        // r is between -1 and 1
    
    int FToV[2 * K][Nv];
    int FToV_T[Nv][2 * K]; // Transposed FToV
    for (int i = 0; i < K; i++) {
        FToV[2 * i][EToV[i][0]], FToV_T[EToV[i][0]][2 * i] = 1;
        FToV[2 * i + 1][EToV[i][1]], FToV_T[EToV[i][1]][2 * i + 1] = 1;
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
} */

void Connect1D(MatrixXd EToV_eigen) {
    int numFaces = 2; // Since this is always 1D
    int totFaces = numFaces * K;
    const int Nv = K + 1;
    int EToV[K][2]; // each row represents element: contains two values from vx for left and right points
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < 2; j++) {
            EToV[i][j] = EToV_eigen(i, j);
        }
    }
    
    int FToV[2 * K][Nv];
    int FToV_T[Nv][2 * K]; // Transposed FToV
    for (int i = 0; i < K; i++) {
        FToV[2 * i][EToV[i][0]], FToV_T[EToV[i][0]][2 * i] = 1;
        FToV[2 * i + 1][EToV[i][1]], FToV_T[EToV[i][1]][2 * i + 1] = 1;
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

    // Basic LGL grid for reference element r
    // VectorXd r(N + 1);
    // int r_length = N + 1;
    r = JacobiGL(0, 0, N);
    cout << "<";
    for (int i = 0; i < N; i++) {
        cout << r(i) << ", ";
    }
    cout << r(N) << ">\n";

    // Reference element matrices
    // MatrixXd V(N + 1, N);
    V = Vandermonde1D(N, r, r_length);
    // MatrixXd invV(N, N + 1);
    invV = V.inverse();

    // MatrixXd Dr(r_length, r_length);
    Dr = Dmatrix1D(r, r_length, N, V);

    // Surface integral terms
    // MatrixXd LIFT(r_length, numFaces * Nfp);
    LIFT = Lift1D(r_length, N, numFaces, Nfp, V);

    // Node Coordinates
    VectorXd va(K);
    VectorXd vb(K);
    va = EToV.col(0).array();
    vb = EToV.col(1).array();
    // MatrixXd x(N + 1, K);
    VectorXd r_intermediary(r_length);
    r_intermediary = r + VectorXd::Ones(r_length);
    x = VectorXd::Ones(N + 1) * VX(va).transpose() + 0.5 * r_intermediary * (VX(vb) - VX(va)).transpose();

    // MatrixXd J(r_length, K);
    J = Jacobian1D(x, Dr, r_length);
    // MatrixXd rx(r_length, K);
    for (int i = 0; i < r_length; i++) {
        for (int j = 0; j < K; j++) {
            rx(i, j) = 1.0 / J(i, j);
        }
    }
    cout << J;
    cout << "\n" << "\n" << rx;

    // Edge Node Masks
    int fmask1_length = 0;
    int fmask2_length = 0;
    for (int i = 0; i < r_length; i++) {
        if ((r(i) + 1.0) < node_tolerance) {
            fmask1_length += 1;
        }
        if ((r(i) - 1.0) < node_tolerance) {
            fmask2_length += 1;
        }
    }
    VectorXd fmask1(fmask1_length);
    VectorXd fmask2(fmask2_length);
    int fm1idx = 0;
    int fm2idx = 0;
    for (int i = 0; i < r_length; i++) {
        if (abs(r(i) + 1.0) < node_tolerance) {
            fmask1(fm1idx) = i;
            fm1idx++;
        }
        if (abs(r(i) - 1.0) < node_tolerance) {
            fmask2(fm2idx) = i;
            fm2idx++;
        }
    } // Something about doing it this way seems really inefficient, but I'm not sure yet how to go about using dynamic array sizes like in python
    MatrixXd fmask(2, fmask1_length); // we're assuming fmask1_length = fmask2_length = 1
    for (int i = 0; i < fmask1_length; i++) {
        fmask(0, i) = fmask1(i);
        fmask(1, i) = fmask2(i);
    }
    MatrixXd Fx(2, K);
    for (int i = 0; i < K; i++) {
        int a = fmask1(0);
        int b = fmask2(0);
        Fx(0, i) = x(a, i);
        Fx(1, i) = x(b, i);
    }

    // Build surface normals and inverse metric at surface
    MatrixXd nx(Nfp * numFaces, K);
    nx = Normals1D(Nfp, numFaces, K);
    MatrixXd Fscale(2, K);
    for (int i = 0; i < K; i++) {
        int a = fmask1(0);
        int b = fmask2(0);
        Fscale(0, i) = 1.0 / J(a, i);
        Fscale(1, i) = 1.0 / J(b, i);
    }

    // Build connectivity matrix
    
}

void AdvecRHS1D(VectorXd u, double time, double a) {
    int alpha = 1;
    MatrixXd du(Nfp * numFaces, K);
    du = MatrixXd::Zero(Nfp * numFaces, K);
    // du(:) = u(vmapM) - u(vmapP)) .* (a * nx(:) ... blah blah

    // Impose boundary condition at x = 0;
    double uin = -sin(a * time);
    // more blah blah

    // Compute right hand sides of the semi-discrete PDE
}

MatrixXd Advec1D(MatrixXd u, double finalTime) {
    MatrixXd u_copy(N + 1, K);
    u_copy = u;
    double time = 0;

    // Runge-Kutta residual storage
    MatrixXd resu(N + 1, K);
    resu = MatrixXd::Zero(N + 1, K);

    // Find minimum space and time step
    double xmin = abs(x(0, 0) - x(1, 0));
    double xmintemp = 0;
    for (int i = 1; i < K; i++) {
        xmintemp = abs(x(0, i) - x(1, i));
        if (xmintemp < xmin) {
            xmin = xmintemp;
        }
    }
    double CFL = 0.75;
    double dt = CFL * xmin / (2 * M_PI);
    dt = dt * 0.5;
    int numSteps = static_cast<int>(ceil(finalTime / dt));
    dt = finalTime / static_cast<double>(numSteps);

    // Advection speed
    double a = 2 * M_PI;

    // Outer time step loop
    for (int tstep = 0; tstep < numSteps; tstep++) {
        for (int INTRK = 0; INTRK < 5; INTRK++) {
            double localTime = time + rk4c(INTRK) * dt;
            // AdvecRHS1D line, but use u_copy
            // resu = rk4a(INTRK) * resu + dt * rhsu;
            // u_copy = u_copy + rk4b(INTRK) * resu;
        }
        // Increment time
        time = time + dt;
    }

    return u_copy;
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
    SelfAdjointEigenSolver<MatrixXd> solver(mat);
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
        PL(1, i) = (static_cast<double>(alpha + beta + 2) * x(i) + static_cast<double>(alpha - beta)) / (2.0 * sqrt(gamma1));
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
        int h1 = 2 * (i + 1) + alpha + beta;
        double anew = 2.0 * sqrt(static_cast<double>((i + 2) * (i + 2 + alpha + beta) * (i + 2 + alpha) * (i + 2 + beta)) / static_cast<double>((h1 + 1) * (h1 + 3))) / (h1 + 2);
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
    } else {
        deltaP = VectorXd::Zero(r_length);
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

MatrixXd Jacobian1D(MatrixXd x, MatrixXd Dr, int r_length) {
   // Used for affine mapping from x to r, rx calculation put in startup
   
   // Calculating Jacobian
   MatrixXd J(r_length, K);
   J = Dr * x;

   return J;
}

MatrixXd Normals1D(int Nfp, int Nfaces, int K) {
    // Compute outward pointing nomrals at element faces

    MatrixXd nx = MatrixXd::Zero(Nfp * Nfaces, K);
    for (int i = 0; i < K; i++) {
        nx(0, i) = -1.0;
        nx(1, i) = 1.0;
    }

    return nx;
}