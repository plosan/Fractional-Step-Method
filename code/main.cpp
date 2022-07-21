#include <iostream>
#include <vector>
#include <cstring>
#include <chrono>
#include "RCGrid.h"
#include "NCMesh.h"
#include "schemes.h"
#include "solver.h"

struct Properties {
    double rho;
    double mu;
};

void printMatrix(const double* a, const int n, const int m, const std::string name);

double* allocateDoubleArray(int n, int m);
void allocateMatrices(double* &u, double* &v, double* &Ru, double* &Rv, const NCMesh m);

// Convection schemes
double schemeUDS(double A, double B, double mf);
double schemeCDS(double A, double B);

// General functions
void computeRuSimplified(double* Ru, const NCMesh m, const double* u, const double* v, const Properties props);
void computeRvSimplified(double* Rv, const NCMesh m, const double* u, const double* v, const Properties props);
void computeRu(double* Ru, const NCMesh m, const double* u, const double* v, const Properties props);
void computeRv(double* Rv, const NCMesh m, const double* u, const double* v, const Properties props);
void computePredictorVelocityU(double* u_pred, const int nx, const int ny, const double* u, const double* Ru, const double* Ru_prev, const Properties props, const double tstep);
void computePredictorVelocityV(double* v_pred, const int nx, const int ny, const double* v, const double* Rv, const double* Rv_prev, const Properties props, const double tstep);
void computeDiscretizationCoefficients(double* A, double* b, const NCMesh m, const double* u_pred, const double* v_pred, const Properties props, const double tstep);
void computeTimeStep(double &tstep, const NCMesh m, const double* u, const double* v, const Properties props);

void allocateOperatorR(const int nx, const int ny, double* &Ru, double* &Rv, double* &Ru_prev, double* &Rv_prev, int& exitCode);

// Lid-driven cavity
namespace lid_driven {
    void setInitialMaps(double* u, double* v, double* p, const NCMesh m, const double u_ref, const double p_ref);
    void setBoundaryPredictorVelocities(double* u_pred, double* v_pred, const NCMesh m, const double u_ref);
    void computeBoundaryDiscretizationCoefficients(double* A, double* b, const NCMesh m);
};




int main(int argc, char* argv[]) {

    double rho = 0.9982;    // Water density at 20 ºC           [kg/m^3]
    double mu = 1.0016e-3;  // Water dynamic viscosity at 20 ºC [Pa s]
    double u_ref = 1;       // X-velocity boundary condition    [m/s]
    double p_ref = 1e5;     // Pressure
    Properties props = {rho, mu};

    const double L = 1;

    const int nx = 10;
    const int ny = 10;
    double tstep = 1e-2;

    const double tol = 1e-12;
    const int maxIt = 500;
    NCMesh m(L, L, 1, nx, ny);
    // m.saveMeshData();
    // m.printMeshData();





    // X-component of velocity
    double* u = (double*) calloc((nx+1)*(ny+2), sizeof(double));
    if(!u) {
        printf("Error: could not allocate enough memory for u\n");
        return false;
    }

    // Y-component of velocity
    double* v = (double*) calloc((nx+2)*(ny+1), sizeof(double));
    if(!v) {
        printf("Error: could not allocate enough memory for v\n");
        return false;
    }

    // Pressure
    double* p = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    if(!p) {
        printf("Error: could not allocate enough memory for p\n");
        return false;
    }

    lid_driven::setInitialMaps(u, v, p, m, u_ref, p_ref);

    double t = 0;   // Current time
    bool steady = false;    // Steady state bool. true = steady state reached, false = steady state not reached

    int exitCode = 0;
    double* Ru;
    double* Rv;
    double* Ru_prev;
    double* Rv_prev;
    allocateOperatorR(nx, ny, Ru, Rv, Ru_prev, Rv_prev, exitCode);
    if(exitCode == -1) {
        printf("Error. The simulation will stop\n");
        return false;
    }

    // printMatrix(Ru, nx+1, ny+2, "Ru");
    // printMatrix(Ru_prev, nx+1, ny+2, "Ru_prev");

    // std::chrono::steady_clock::time_point begin, end;


    computeRu(Ru, m, u, v, props);
    computeRuSimplified(Ru_prev, m, u, v, props);
    std::memcpy(Ru_prev, Ru, (nx+1)*(ny+2)*sizeof(double));

    computeRvSimplified(Rv, m, u, v, props);
    computeRvSimplified(Rv_prev, m, u, v, props);
    std::memcpy(Rv_prev, Rv, (nx+2)*(ny+1)*sizeof(double));

    double* u_pred = (double*) calloc((nx+1)*(ny+2), sizeof(double));
    if(!u_pred) {
        printf("Error: could not allocate enough memory for u_pred\n");
        return false;
    }

    double* v_pred = (double*) calloc((nx+2)*(ny+1), sizeof(double));
    if(!v_pred) {
        printf("Error: could not allocate enough memory for u_pred\n");
        return false;
    }

    computePredictorVelocityU(u_pred, nx, ny, u, Ru, Ru_prev, props, tstep);
    computePredictorVelocityV(v_pred, nx, ny, v, Rv, Rv_prev, props, tstep);

    lid_driven::setBoundaryPredictorVelocities(u_pred, v_pred, m, u_ref);

    double* A = (double*) calloc(5*(nx+2)*(ny+2), sizeof(double));
    if(!A) {
        printf("Error: could not allocate enough memory for coefficients matrix\n");
        return false;
    }

    double* b = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    if(!b) {
        printf("Error: could not allocate enough memory for independent terms vector\n");
        return false;
    }


    computeDiscretizationCoefficients(A, b, m, u_pred, v_pred, props, tstep);
    lid_driven::computeBoundaryDiscretizationCoefficients(A, b, m);

    for(int j = 0; j < ny+2; j++) {
        for(int i = 0; i < nx+2; i++) {
            int k = j * (nx + 2) + i;
            printf("(%3d,%3d)%10d", i, j, k);
            for(int i = 0; i < 5; i++)
                printf("%15.2f", A[5*k+i]);
            printf("%15.2f\n", b[k]);
        }
    }

    // for(int k = 0; k < (nx+2)*(ny+2); k++) {
    //     printf("%5d", k);
    //     for(int i = 0; i < 5; i++)
    //         printf("%15.2f", A[5*k+i]);
    //     printf("%15.2f\n", b[k]);
    // }

    printMatrix(p, nx+2, ny+2, "p0");
    int exitCodeGS = solveSystemGS(nx+2, ny+2, tol, maxIt, A, b, p);
    printMatrix(p, nx+2, ny+2, "p");

    printf("exitCode: %d\n", exitCodeGS);

}



void printMatrix(const double* a, const int n, const int m, const std::string name) {

    // Print name
    printf("\n%s = \n", name.c_str());

    // Table header
    printf("%10s", "");
    for(int i = 0; i < n; i++)
        printf("%10d", i);
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

void computeRuSimplified(double* Ru, const NCMesh m, const double* u, const double* v, const Properties props) {


    int nx = m.getNX();
    int ny = m.getNY();

    for(int j = 1; j < ny+1; j++) {

        for(int i = 1; i < nx; i++) {

            // Nodal velocities
            int node = j*(nx+1) + i;
            double uP = u[node];
            double uW = u[node-1];
            double uE = u[node+1];
            double uS = u[node-(nx+1)];
            double uN = u[node+(nx+1)];

            // Faces velocities
            double uw = schemeCDS(uP, uW);
            double ue = schemeCDS(uP, uE);
            double us = (j > 1 ? schemeCDS(uP, uS) : uS);
            double un = (j < ny ? schemeCDS(uP, uN) : uN);

            // Areas
            double Ax = m.atSurfX(j);
            double Ay = m.atSurfY_StaggX(i);
            double Ay_left = m.atSemiSurfY(i,1);
            double Ay_right = m.atSemiSurfY(i+1,0);

            // Mass flows
            double mw = 0.5 * props.rho * (uP + uW) * Ax;
            double me = 0.5 * props.rho * (uP + uE) * Ax;
            double mn = props.rho * (v[j*(nx+2)+i] * Ay_left + v[j*(nx+2)+i+1] * Ay_right);
            double ms = props.rho * (v[(j-1)*(nx+2)+i] * Ay_left + v[(j-1)*(nx+2)+i+1] * Ay_right);

            // Operator R(u)
            double integral1 = -(me * ue - mw * uw + mn * un - ms * us);
            double integral2 = Ax * (uE - uP) / m.atDistFaceX(i) - Ax * (uP - uW) / m.atDistFaceX(i-1);
            integral2 += Ay * (uN - uP) / m.atDistY(j) - Ay * (uP - uS) / m.atDistY(j-1);
            integral2 *= props.mu;
            Ru[j*(nx+1)+i] = (integral1 + integral2) / m.atVolStaggX(i,j);


        }

    }

}

void computeRu(double* Ru, const NCMesh m, const double* u, const double* v, const Properties props) {

    int nx = m.getNX();
    int ny = m.getNY();

    // To save the mass flows of the south face
    // The mass flow through the north face of node (i,j) equals the mass flow through the south face of node (i,j+1)
    // When the first row is computed (j = 1), my is filled with the north face mass flows. In subsequent rows (j > 1), the mass flow through the
    // south face of node (i,j) is taken from my[i] and then my[i] is replaced by the mass flow through the north face of node (i,j)
    double* my = (double*) calloc(nx+1, sizeof(double));
    if(!my) {
        printf("Error: could not allocate enough memory for variable my\n");
        return;
    }

    // First row (j = 1)
    int j = 1;

    // Properties that can be updated instead of computing them twice
    double uW = u[j*(nx+1)];        // West node velocity
    double uP = u[j*(nx+1)+1];      // Current node velocity
    double uw = schemeCDS(uW, uP);  // Velocity west face
    double mw = 0.5 * props.rho * (uW + uP) * m.atSurfX(j); // Mass flow west face

    for(int i = 1; i < nx; i++) {
        // Nodal velocities
        double uE = u[j*(nx+1)+i+1];    // East node velocity
        double uS = u[(j-1)*(nx+1)+i];  // South node velocity
        double uN = u[(j+1)*(nx+1)+i];  // North node velocity

        // Mass flows
        double me = 0.5 * props.rho * (uP + uE) * m.atSurfX(j);                                                         // Mass flow east face
        double mn = props.rho * (v[j*(nx+2)+i] * m.atSemiSurfY(i,1) + v[j*(nx+2)+i+1] * m.atSemiSurfY(i+1,0));          // Mass flow north face
        double ms = props.rho * (v[(j-1)*(nx+2)+i] * m.atSemiSurfY(i,1) + v[(j-1)*(nx+2)+i+1] * m.atSemiSurfY(i+1,0));  // Mass flow south face
        my[i] = mn;

        // Face velocities
        double ue = schemeCDS(uP, uE);
        double us = uS;                 // Equal to uS because it is the first row
        double un = schemeCDS(uP, uN);  // This may give problems

        // Operator R(u)
        double integral1 = -(me * ue - mw * uw + mn * un - ms * us);
        double integral2 = m.atSurfX(j) * (uE - uP) / m.atDistFaceX(i) - m.atSurfX(j) * (uP - uW) / m.atDistFaceX(i-1);
        integral2 += m.atSurfY_StaggX(i) * (uN - uP) / m.atDistY(j) - m.atSurfY_StaggX(i) * (uP - uS) / m.atDistY(j-1);
        integral2 *= props.mu;
        Ru[j*(nx+1)+i] = (integral1 + integral2) / m.atVolStaggX(i,j);

        // Next node
        uW = uP;
        uP = uE;
        uw = ue;
        mw = me;
    }

    // In between rows
    for(j = 2; j < ny; j++) {
        // Properties that can be updated instead of computing them twice
        uW = u[j*(nx+1)];                                   // West node velocity
        uP = u[j*(nx+1)+1];                                 // Current node velocity
        uw = schemeCDS(uW, uP);                             // Velocity west face
        mw = 0.5 * props.rho * (uW + uP) * m.atSurfX(j);    // Mass flow west face
        // Go through all nodes in the j-th row
        for(int i = 1; i < nx; i++) {
            // Nodal velocities
            double uE = u[j*(nx+1)+i+1];
            double uS = u[(j-1)*(nx+1)+i];
            double uN = u[(j+1)*(nx+1)+i];

            // Mass flows
            double me = 0.5 * props.rho * (uP + uE) * m.atSurfX(j);                                                 // East face mass flow
            double ms = my[i];                                                                                      // South face mass flow (previously computed as mass flow through north face)
            double mn = props.rho * (v[j*(nx+2)+i] * m.atSemiSurfY(i,1) + v[j*(nx+2)+i+1] * m.atSemiSurfY(i+1,0));  // North face mass flow
            my[i] = mn;                                                                                             // Save north face mass flow (to later on be used as mass flow through south face)

            // Face velocities
            double ue = schemeCDS(uP, uE);
            double us = schemeCDS(uP, uS);  // This may give problems. It can be optimised by computing it in the previous row
            double un = schemeCDS(uP, uN);  // This may give problems

            // Operator R(u)
            double integral1 = -(me * ue - mw * uw + mn * un - ms * us);
            double integral2 = m.atSurfX(j) * (uE - uP) / m.atDistFaceX(i) - m.atSurfX(j) * (uP - uW) / m.atDistFaceX(i-1);
            integral2 += m.atSurfY_StaggX(i) * (uN - uP) / m.atDistY(j) - m.atSurfY_StaggX(i) * (uP - uS) / m.atDistY(j-1);
            integral2 *= props.mu;
            Ru[j*(nx+1)+i] = (integral1 + integral2) / m.atVolStaggX(i,j);

            // Next node
            uW = uP;
            uP = uE;
            uw = ue;
            mw = me;
        }
    }

    // Last row (j = ny)
    j = ny;

    uW = u[j*(nx+1)];
    uP = u[j*(nx+1)+1];
    uw = schemeCDS(uW, uP);
    mw = 0.5 * props.rho * (uW + uP) * m.atSurfX(j);

    for(int i = 1; i < nx; i++) {
        // Nodal velocities
        double uE = u[j*(nx+1)+i+1];   // East node velocity
        double uS = u[(j-1)*(nx+1)+i]; // South node velocity
        double uN = u[(j+1)*(nx+1)+i]; // North node velocity

        // Mass flows
        double me = 0.5 * props.rho * (uP + uE) * m.atSurfX(j);
        double ms = my[i];
        double mn = props.rho * (v[(j+1)*(nx+2)+i] * m.atSemiSurfY(i,1) + v[(j+1)*(nx+2)+i+1] * m.atSemiSurfY(i+1,0));

        // Face velocities
        double ue = schemeCDS(uP, uE);
        double us = schemeCDS(uP, uS);  // This may give problems
        double un = uN;                 // Equal to uN because it is the last row

        // Operator Ru
        double integral1 = -(me * ue - mw * uw + mn * un - ms * us);
        double integral2 = m.atSurfX(j) * (uE - uP) / m.atDistFaceX(i) - m.atSurfX(j) * (uP - uW) / m.atDistFaceX(i-1);
        integral2 += m.atSurfY_StaggX(i) * (uN - uP) / m.atDistY(j) - m.atSurfY_StaggX(i) * (uP - uS) / m.atDistY(j-1);
        integral2 *= props.mu;
        Ru[j*(nx+1)+i] = (integral1 + integral2) / m.atVolStaggX(i,j);

        // Next node
        uW = uP;
        uP = uE;
        uw = ue;
        mw = me;
    }
}

void computeRvSimplified(double* Rv, const NCMesh m, const double* u, const double* v, const Properties props) {

    int nx = m.getNX();
    int ny = m.getNY();

    for(int i = 1; i < nx+1; i++) {

        for(int j = 1; j < ny; j++) {

            int node = j * (nx + 2) + i;
            // Nodal velocities
            double vP = v[node];
            double vW = v[node-1];
            double vE = v[node+1];
            double vS = v[node-(nx+2)];
            double vN = v[node+(nx+2)];

            // Face velocities
            double vw = (i > 1 ? schemeCDS(vP, vW) : vW);
            double ve = (i < nx ? schemeCDS(vP, vE) : vE);
            double vs = schemeCDS(vP, vS);
            double vn = schemeCDS(vP, vN);

            // Areas
            double Ax = m.atSurfX_StaggY(j);
            double Ax_up = m.atSemiSurfX(j+1,0);
            double Ax_down = m.atSemiSurfX(j,1);
            double Ay = m.atSurfY(i);

            // Mass flows
            double mw = props.rho * (u[j*(nx+1)+i-1] * Ax_down + u[(j+1)*(nx+1)+i-1] * Ax_up);
            double me = props.rho * (u[j*(nx+1)+i] * Ax_down + u[(j+1)*(nx+1)+i+1] * Ax_up);
            double ms = 0.5 * props.rho * (vP + vS) * Ay;
            double mn = 0.5 * props.rho * (vP + vN) * Ay;

            // Operator R(v)
            double integral1 = -(me * ve - mw * vw + mn * vn - ms * vs);
            double integral2 = Ax * (vE - vP) / m.atDistX(i) - Ax * (vP - vW) / m.atDistX(i-1);
            integral2 += Ay * (vN - vP) / m.atDistFaceY(j) - Ay * (vP - vS) / m.atDistFaceY(j-1);
            integral2 *= props.mu;
            Rv[j*(nx+2)+i] = (integral1 + integral2) / m.atVolStaggY(i,j);

        }

}

}

void computeRv(double* Rv, const NCMesh m, const double* u, const double* v, const Properties props) {


    int nx = m.getNX();
    int ny = m.getNY();

    double* mx = (double*) calloc(ny+1, sizeof(double));
    if(!mx) {
        printf("Error: could not allocate enough memory for variable mx\n");
        return;
    }

    // Column i = 1
    int i = 1;
    int j = 1;
    double vS = v[(j-1)*(nx+2)+i];                          // South node y-velocity
    double vP = v[j*(nx+2)+i];                              // Current node y-velocity
    double vs = schemeCDS(vP, vS);                          // South face y-velocity
    double ms = 0.5 * props.rho * (vP + vS) * m.atSurfY(i); // South face mass flow

    for(j = 1; j < ny; j++) {

        // Nodal velocities
        double vW = v[j*(nx+2)+i-1];    // West node y-velocity
        double vE = v[j*(nx+2)+i+1];    // East node y-velocity
        double vN = v[(j+1)*(nx+2)+i];  // North node y-velocity

        // Mass flows
        double mw = props.rho * (u[j*(nx+1)+i-1] * m.atSemiSurfX(j,1) + u[(j+1)*(nx+1)+i-1] * m.atSemiSurfX(j+1,0));    // West face mass flow
        double me = props.rho * (u[j*(nx+1)+i] * m.atSemiSurfX(j,1) + u[(j+1)*(nx+1)+i] * m.atSemiSurfX(j+1,0));        // East face mass flow
        double mn = 0.5 * props.rho * (vP + vN) * m.atSurfY(i);                                                         // North face mass flow
        mx[j] = me;

        // Face velocities
        double vw = vW;                 // West face y-velocity
        double ve = schemeCDS(vP, vE);  // East face y-velocity. This may give problems
        double vn = schemeCDS(vP, vN);  // North face y-velocity

        // Operator R(v)
        double integral1 = -(me * ve - mw * vw + mn * vn - ms * vs);
        double integral2 = m.atSurfX_StaggY(j) * (vE - vP) / m.atDistX(i) - m.atSurfX_StaggY(j) * (vP - vW) / m.atDistX(i-1);
        integral2 += m.atSurfY(i) * (vN - vP) / m.atDistFaceY(j) - m.atSurfY(i) * (vP - vS) / m.atDistFaceY(j-1);
        integral2 *= props.mu;
        Rv[j*(nx+2)+i] = (integral1 + integral2) / m.atVolStaggY(i,j);

        // Next node
        vS = vP;
        vP = vE;
        vs = vn;
        ms = mn;

    }

    // Inner columns (2 <= i < nx)
    for(i = 2; i < nx; i++) {

        vS = v[(j-1)*(nx+2)+i];                             // South node y-velocity
        vP = v[j*(nx+2)+i];                                 // Current node y-velocity
        vs = schemeCDS(vP, vS);                             // South face y-velocity
        ms = 0.5 * props.rho * (vP + vS) * m.atSurfY(i);    // South face mass flow

        for(int j = 1; j < ny; j++) {

            // Nodal velocities
            double vW = v[j*(nx+2)+i-1];
            double vE = v[j*(nx+2)+i+1];
            double vN = v[(j+1)*(nx+2)+i];

            // Mass flows
            double mw = mx[j];                                                                                          // West face mass flow
            double me = props.rho * (u[j*(nx+1)+i] * m.atSemiSurfX(j,1) + u[(j+1)*(nx+1)+i] * m.atSemiSurfX(j+1,0));    // West face mass flow
            double mn = 0.5 * props.rho * (vP + vN) * m.atSurfY(i);                                                     // North face mass flow
            mx[j] = me;

            // Face velocities
            double vw = schemeCDS(vP, vW);  // West face y-velocity. This may give problems
            double ve = schemeCDS(vP, vE);  // East face y-velocity. This may give problems
            double vn = schemeCDS(vP, vN);  // North face y-velocity. This may give problems

            // Operator R(v)
            double integral1 = -(me * ve - mw * vw + mn * vn - ms * vs);
            double integral2 = m.atSurfX_StaggY(j) * (vE - vP) / m.atDistX(i) - m.atSurfX_StaggY(j) * (vP - vW) / m.atDistX(i-1);
            integral2 += m.atSurfY(i) * (vN - vP) / m.atDistFaceY(j) - m.atSurfY(i) * (vP - vS) / m.atDistFaceY(j-1);
            integral2 *= props.mu;
            Rv[j*(nx+2)+i] = (integral1 + integral2) / m.atVolStaggY(i,j);

            // Next node
            vS = vP;
            vP = vN;
            vs = vn;
            ms = mn;

        }


    }

    // Last column (i = nx)
    i = nx;

    for(j = 1; j < ny; j++) {

        // Nodal velocities
        double vW = v[j*(nx+2)+i-1];    // West node y-velocity
        double vE = v[j*(nx+2)+i+1];    // East node y-velocity
        double vN = v[(j+1)*(nx+2)+i];  // North node y-velocity

        // Mass flows
        double mw = mx[j];
        double me = props.rho * (u[j*(nx+1)+i] * m.atSemiSurfX(j,1) + u[(j+1)*(nx+1)+i] * m.atSemiSurfX(j+1,0));
        double mn = 0.5 * props.rho * (vP + vN) * m.atSurfY(i);

        // Face velocities
        double vw = schemeCDS(vP, vW);  // West face y-velocity. This may give problems
        double ve = schemeCDS(vP, vE);  // East face y-velocity. This may give problems
        double vn = schemeCDS(vP, vN);  // North face y-velocity

        // Operator R(v)
        double integral1 = -(me * ve - mw * vw + mn * vn - ms * vs);
        double integral2 = m.atSurfX_StaggY(j) * (vE - vP) / m.atDistX(i) - m.atSurfX_StaggY(j) * (vP - vW) / m.atDistX(i-1);
        integral2 += m.atSurfY(i) * (vN - vP) / m.atDistFaceY(j) - m.atSurfY(i) * (vP - vS) / m.atDistFaceY(j-1);
        integral2 *= props.mu;
        Rv[j*(nx+2)+i] = (integral1 + integral2) / m.atVolStaggY(i,j);

        // Next node
        vS = vP;
        vP = vN;
        vs = vn;
        ms = mn;
    }



}

void computePredictorVelocityU(double* u_pred, const int nx, const int ny, const double* u, const double* Ru, const double* Ru_prev, const Properties props, const double tstep) {

    // int nx = m.getNX();
    // int ny = m.getNY();

    for(int i = 1; i < nx; i++) {
        for(int j = 1; j < ny+1; j++) {
            int k = j * (nx + 1) + i;
            u_pred[k] = u[k] + tstep / props.rho * (1.5 * Ru[k] - 0.5 * Ru_prev[k]);
        }
    }

}


void computePredictorVelocityV(double* v_pred, const int nx, const int ny, const double* v, const double* Rv, const double* Rv_prev, const Properties props, const double tstep) {

    // int nx = m.getNX();
    // int ny = m.getNY();

    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny; j++) {
            int k = j*(nx+2) + i;
            v_pred[k] = v[k] + tstep / props.rho * (1.5 * Rv[k] - 0.5 * Rv_prev[k]);
        }
    }

}

void computeDiscretizationCoefficients(double* A, double* b, const NCMesh m, const double* u_pred, const double* v_pred, const Properties props, const double tstep) {

    int nx = m.getNX();
    int ny = m.getNY();

    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny+1; j++) {
            int node = j * (nx + 2) + i;

            double Ax = m.atSurfX(j);
            double Ay = m.atSurfY(i);
            double dPE = m.atDistX(i);
            double dPW = m.atDistX(i-1);
            double dPN = m.atDistY(j);
            double dPS = m.atDistY(j-1);

            // South node
            A[5*node] = Ay / dPS;

            // West node
            A[5*node+1] = Ax / dPW;

            // East node
            A[5*node+2] = Ax / dPE;

            // North node
            A[5*node+3] = Ay / dPN;

            // Current node
            A[5*node+4] = A[5*node] + A[5*node+1] + A[5*node+2] + A[5*node+3];

            // Independent term
            int node_x = j * (nx + 1) + i;
            int node_y = j * (nx + 2) + i;
            b[node] = u_pred[node_x] * Ax - u_pred[node_x-1] * Ax;          // Independent term: predictor x-velocity terms
            b[node] += v_pred[node_y] * Ay - v_pred[node_y-(nx+2)] * Ay;    // Independent term: predictor y-velocity terms
            b[node] *= (-props.rho) / tstep;

        }
    }

}

void computeTimeStep(double &tstep, const NCMesh m, const double* u, const double* v, const Properties props) {

    int nx = m.getNX();
    int ny = m.getNY();

    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny+1; j++) {

        }
    }

}

void allocateOperatorR(const int nx, const int ny, double* &Ru, double* &Rv, double* &Ru_prev, double* &Rv_prev, int &exitCode) {

    exitCode = -1;

    // Operator Ru at time n
    Ru = (double*) calloc((nx+1)*(ny+2), sizeof(double));
    if(!Ru) {
        printf("Error: could not allocate enough memory for Ru\n");
        return;
    }

    // Operator Rv at time n
    Rv = (double*) calloc((nx+2)*(ny+1), sizeof(double));
    if(!Rv) {
        printf("Error: could not allocate enough memory for Rv\n");
        return;
    }

    // Operator Ru at time n-1
    Ru_prev = (double*) calloc((nx+1)*(ny+2), sizeof(double));
    if(!Ru_prev) {
        printf("Error: could not allocate enough memory for Ru\n");
        return;
    }

    // Operator Rv at time n-1
    Rv_prev = (double*) calloc((nx+2)*(ny+1), sizeof(double));
    if(!Rv) {
        printf("Error: could not allocate enough memory for Rv_prev\n");
        return;
    }

    exitCode = 1;   // All arrays were allocated successfully
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


// Lid-driven cavity

void lid_driven::setInitialMaps(double* u, double* v, double* p, const NCMesh m, const double u_ref, const double p_ref) {

    int nx = m.getNX(); // Number of control volumes along x axis
    int ny = m.getNY(); // Number of control volumes along y axis

    // Arrays u, v and p allocated using calloc, all bits set to zero by calloc

    // x-velocity (u) initial map
    for(int i = 0; i < nx+1; i++)
        u[(ny+1)*(nx+1)+i] = u_ref;

    // pressure (p) initial map
    for(int k = 0; k < (nx+2)*(ny+2); k++)
        p[k] = p_ref;
}


void lid_driven::setBoundaryPredictorVelocities(double* u_pred, double* v_pred, const NCMesh m, const double u_ref) {

    int nx = m.getNX();
    int ny = m.getNY();

    // PREDICTOR VELOCITY X
    // Lower and upper boundaries
    for(int i = 0; i < nx+1; i++) {
        u_pred[i] = 0;                      // Lower boundary
        u_pred[(ny+1)*(nx+1)+i] = u_ref;    // Upper boundary
    }

    // Left and right boundaries
    for(int j = 1; j < ny+2; j++) {
        u_pred[j*(nx+1)] = 0;           // Left boundary
        u_pred[j*(nx+1)+nx] = 0;        // Right boundary
    }

    // PREDICTOR VELOCITY Y
    // Lower and upper boundaries
    for(int i = 0; i < nx+2; i++) {
        v_pred[i] = 0;                  // Lower boundary
        v_pred[ny*(nx+2)+i] = 0;        // Upper boundary
    }

    // Left and right boundaries
    for(int j = 1; j < ny+1; j++) {
        v_pred[j*(nx+2)] = 0;           // Left boundary
        v_pred[j*(nx+2)+nx+1] = 0;      // Right boundary
    }

}

void lid_driven::computeBoundaryDiscretizationCoefficients(double* A, double* b, const NCMesh m) {

    int nx = m.getNX();
    int ny = m.getNY();

    // lOWER AND UPPER BOUNDARIES (including corner nodes, i.e. singular nodes)
    for(int i = 0; i < nx+2; i++) {
        // Lower boundary
        int node = i;       // Node identifier
        A[5*node] = 0;      // South node
        A[5*node+1] = 0;    // West node
        A[5*node+2] = 0;    // East node
        A[5*node+3] = 1;    // North node
        A[5*node+4] = 1;    // Central node
        b[node] = 0;        // Independent term
        // Upper boundary
        node = (ny + 1) * (nx + 2) + i;
        A[5*node] = 1;      // South node
        A[5*node+1] = 0;    // West node
        A[5*node+2] = 0;    // East node
        A[5*node+3] = 0;    // North node
        A[5*node+4] = 1;    // Central node
        b[node] = 0;        // Independent term
    }

    // LEFT AND RIGHT BOUNDARIES
    for(int j = 1; j < ny+1; j++) {
        // Left boundary
        int node = j * (nx + 2);
        A[5*node] = 0;      // South node
        A[5*node+1] = 0;    // West node
        A[5*node+2] = 1;    // East node
        A[5*node+3] = 0;    // North node
        A[5*node+4] = 1;    // Central node
        b[node] = 0;        // Independent term
        // Right boundary
        node = j * (nx + 2) + (nx + 1);
        A[5*node] = 0;      // South node
        A[5*node+1] = 1;    // West node
        A[5*node+2] = 0;    // East node
        A[5*node+3] = 0;    // North node
        A[5*node+4] = 1;    // Central node
        b[node] = 0;        // Independent term
    }
}
