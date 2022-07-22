#include <iostream>
#include <vector>
#include <cstring>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "NCMesh.h"
// #include "schemes.h"
#include "solver.h"
#include "fsm.h"



void printMatrix(const double* a, const int n, const int m, const std::string name);
void printLinearSystem(const double* A, const double* b, const int nx, const int ny);

double* allocateDoubleArray(int n, int m);

// Convection schemes
double schemeUDS(double A, double B, double mf);
double schemeCDS(double A, double B);

int main(int argc, char* argv[]) {

    double rho = 1;    // Density
    double mu = 1e-4;  // Dynamic viscosity
    double u_ref = 1e-2;       // X-velocity boundary condition    [m/s]
    double p_ref = 1e5;     // Pressure
    Properties props = {rho, mu};

    const double L = 1;

    const int nx = 129;
    const int ny = 129;
    double tstep = 1e-2;

    const double tol = 1e-12;       // Linear system solver tolerance
    const int maxIt = 1e5;         // Linear system max iterations
    const double sstol = 1e-7;      // Tolerance for steady state check
    NCMesh m(L, L, 1, nx, ny);


    double t = 0;           // Current time
    int it = 0;             // Iteration counter
    bool steady = false;    // Steady state bool. true = steady state reached, false = steady state not reached

    double* u;
    double* v;
    double* p;
    allocateFluidVariables(nx, ny, u, v, p);
    lid_driven::setInitialMaps(u, v, p, m, u_ref, p_ref);

    double* Ru;
    double* Rv;
    double* Ru_prev;
    double* Rv_prev;
    allocateOperatorR(nx, ny, Ru, Rv, Ru_prev, Rv_prev);

    computeRu(Ru_prev, m, u, v, props);
    computeRv(Rv_prev, m, u, v, props);

    double* u_pred;
    double* v_pred;
    allocatePredictorVelocities(nx, ny, u_pred, v_pred);

    double* A;
    double* b;
    allocateLinearSystemVariables(nx, ny, A, b);

    std::chrono::steady_clock::time_point begin, end;

    while(!steady) {

        begin = std::chrono::steady_clock::now();

        // Compute operator R(u) and R(v)
        computeRu(Ru, m, u, v, props);
        computeRv(Rv, m, u, v, props);

        // Compute predictor velocities
        computePredictorVelocityU(u_pred, nx, ny, u, Ru, Ru_prev, props, tstep);
        computePredictorVelocityV(v_pred, nx, ny, v, Rv, Rv_prev, props, tstep);
        lid_driven::setBoundaryPredictorVelocities(u_pred, v_pred, nx, ny, u_ref);

        // Compute discretization coefficients
        computeDiscretizationCoefficients(A, b, m, u_pred, v_pred, props, tstep);
        lid_driven::computeBoundaryDiscretizationCoefficients(A, b, m);

        // Solve linear system
        int exitCodeGS = solveSystemGS(nx+2, ny+2, tol, maxIt, A, b, p);
        if(exitCodeGS == -1) {
            printf("Error solviusaodhsadisadsa\n");
            return false;
        }

        double maxDiff = 0;
        computeVelocity(u, v, m, u_pred, v_pred, p, props, tstep, maxDiff);

        t += tstep;
        it++;

        end = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1e6;


        printf("%6d %5s %10.5f %5s %10.5f %5s %10.5e %5s %10d %5s %.3f\n", it, "", t, "", tstep, "", maxDiff, "", exitCodeGS, "", elapsed);

        if(it % 100 == 0) {

            double* u_col = (double*) calloc((nx+2)*(ny+2), sizeof(double));
            double* v_col = (double*) calloc((nx+2)*(ny+2), sizeof(double));
            if(!u_col) {
                printf("Error: could not allocate enough memory for u_col\n");
                return false;
            }
            if(!v_col) {
                printf("Error: could not allocate enough memory for v_col\n");
                return false;
            }
            computeVelocityCollocatedMesh(u_col, v_col, nx, ny, u, v);

            // printMatrix(u_col, nx+2, ny+2, "u_col");
            // printMatrix(v_col, nx+2, ny+2, "v_col");

            printVelocityToFile(m, u_col, v_col, "sim/vel.txt", 5);
            printVelocityUToFile(m, u_col, "sim/u.txt", 5);
            printVelocityVToFile(m, v_col, "sim/v.txt", 5);
            printPressureToFile(m, p, "sim/p.txt", 5);

            free(u_col);
            free(v_col);

        }

        if(maxDiff < sstol)
            steady = true;
        else {
            // Compute next time step
            computeTimeStep(tstep, m, u, v, props, tol);

            // Update operator R
            std::memcpy(Ru_prev, Ru, (nx+1)*(ny+2)*sizeof(double));
            std::memcpy(Rv_prev, Rv, (nx+2)*(ny+1)*sizeof(double));
        }

        // printf("here\n");
    }



    // printMatrix(p, nx+2, ny+2, "p");
    // printMatrix(u, nx+1, ny+2, "u");
    // printMatrix(v, nx+2, ny+1, "v");



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










double* allocateDoubleArray(int n, int m) {
    /*
    allocateDoubleArray: allocates an array of doubles sequentially
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - n     Row count
        - m     Column count
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - double*   Double array of size n*m allocated sequentially
    */
    double* a = (double*) calloc(n*m, sizeof(double));
    return a;
}

double schemeUDS(double A, double B, double mf) {
    /*
    schemeUDS: computes the value of the variable according to the Upwind-Difference Scheme (UDS)
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - A     Value of the variable downstream (with respect to the face where the convective property is being computed)
        - B     Value of the variable upstream (with respect to the face where the convective property is being computed)
        - mf    Mass flow through through the face where the convective property is being computed
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - double    Convection variable
    */
    if(mf > 0)
        return A;
    if(mf < 0)
        return B;
    return schemeCDS(A, B);
}

double schemeCDS(double A, double B) {
    return 0.5*(A + B);
}
