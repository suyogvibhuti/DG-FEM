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