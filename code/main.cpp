#include <iostream>

// #include <algorithm>
// #include <chrono>
// #include <cmath>
// #include <cstring>
// #include <fstream>
// #include <iomanip>
// #include <vector>

#include "fsm.h"



void printMatrix(const double* a, const int n, const int m, const std::string name);
void printLinearSystem(const double* A, const double* b, const int nx, const int ny);

int main(int argc, char* argv[]) {

    double rho = 1;    // Density
    double mu = 1e-4;  // Dynamic viscosity
    double u_ref = 1e-1;       // X-velocity boundary condition    [m/s]
    double p_ref = 1e5;     // Pressure

    const double L = 1;

    const int nx = 129;
    const int ny = 129;
    double tstep = 1e-3;

    const double tol = 1e-12;       // Linear system solver tolerance
    const int maxIt = 1e6;         // Linear system max iterations
    const double sstol = 1e-4;      // Tolerance for steady state check


    lid_driven::mainLoop(rho, mu, u_ref, p_ref, L, nx, ny, tstep, maxIt, tol, sstol);

}



void printMatrix(const double* a, const int n, const int m, const std::string name) {

    // Print name
    printf("\n%s = \n", name.c_str());

    // Table header
    printf("%10s", "");
    for(int i = 0; i < n; i++)
        printf("%12d", i);
    printf("\n");

    // Table content
    for(int j = m-1; j >= 0; j--) {
        printf("%5d%5s", j, "");
        for(int i = 0; i < n; i++)
            printf("%12.2f", a[j*n+i]);
        printf("\n");
    }
    printf("\n");

}

void printLinearSystem(const double* A, const double* b, const int nx, const int ny) {
    for(int j = 0; j < ny+2; j++) {
        for(int i = 0; i < nx+2; i++) {
            int k = j * (nx + 2) + i;
            printf("(%3d,%3d)%10d", i, j, k);
            for(int i = 0; i < 5; i++)
                printf("%15.2f", A[5*k+i]);
            printf("%15.2f\n", b[k]);
        }
    }
}
