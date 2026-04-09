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
    // Running dgfem function, dummy value for fluid velocity
    dgfem(2.0, 32.0);

    // Exiting program, normal operation code
    return 0;
}

void dgfem(double fluidVelocity, double length) {
    cout << "hello!";
}