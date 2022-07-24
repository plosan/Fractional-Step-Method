#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cfloat>
#include <cstring>
#include <chrono>

#include "fsm.h"
#include "schemes.h"
#include "solver.h"


void allocateFluidVariables(const int nx, const int ny, double* &u, double* &v, double* &p) {
    /*
    allocateFluidVariables: allocates the arrays for x-velocity, y-velocity and pressure
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx    Control volumes in x direction          [const int]
        - ny    Control volumes in y direction          [const int]
        - u     X-velocity array. Size: (nx+1)*(ny+2)   [double* &]
        - v     Y-velocity array. Size: (nx+2)*(ny+1)   [double* &]
        - p     Pressure. Size: (nx+2)*(ny+2)           [double* v]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - none
    */
    // X-component of velocity
    u = (double*) calloc((nx+1)*(ny+2), sizeof(double));
    if(!u) {
        printf("Error: could not allocate enough memory for u\n");
        return;
    }
    // Y-component of velocity
    v = (double*) calloc((nx+2)*(ny+1), sizeof(double));
    if(!v) {
        printf("Error: could not allocate enough memory for v\n");
        return;
    }
    // Pressure
    p = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    if(!p) {
        printf("Error: could not allocate enough memory for p\n");
        return;
    }
}

void allocateOperatorR(const int nx, const int ny, double* &Ru, double* &Rv, double* &Ru_prev, double* &Rv_prev) {
    /*
    allocateOperatorR: allocates the arrays for the operator R, namely Ru, Rv, Ru_prev and Rv_prev
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx        Control volumes in x direction                                          [const int]
        - ny        Control volumes in y direction                                          [const int]
        - Ru        X-component of operator R at the current time. Size: (nx+1)*(ny+2)      [double* &]
        - Rv        Y-component of operator R at the current time. Size: (nx+2)*(ny+1)      [double* &]
        - Ru_prev   X-component of operator R at the previous time. Size: (nx+1)*(ny+2)     [double* &]
        - Rv_prev   Y-component of operator R at the previous time. Size: (nx+2)*(ny+1)     [double* &]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - none
    */
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
}


void allocatePredictorVelocities(const int nx, const int ny, double* &u_pred, double* &v_pred) {
    // X-velocity predictor
    u_pred = (double*) calloc((nx+1)*(ny+2), sizeof(double));
    if(!u_pred) {
        printf("Error: could not allocate enough memory for u_pred\n");
        return;
    }
    // Y-velocity predictor
    v_pred = (double*) calloc((nx+2)*(ny+1), sizeof(double));
    if(!v_pred) {
        printf("Error: could not allocate enough memory for u_pred\n");
        return;
    }
}


void allocateLinearSystemVariables(const int nx, const int ny, double* &A, double* &b) {
    // Linear system matrix
    A = (double*) calloc(5*(nx+2)*(ny+2), sizeof(double));
    if(!A) {
        printf("Error: could not allocate enough memory for coefficients matrix\n");
        return;
    }
    // Linear system vector
    b = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    if(!b) {
        printf("Error: could not allocate enough memory for independent terms vector\n");
        return;
    }
}


void computeMassFlowsStaggX(double* mx, double* my, const NCMesh m, const double* u, const double* v, const Properties props) {

    int nx = m.getNX();
    int ny = m.getNY();

    mx = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    for(int j = 1; j < ny+1; j++) {
        double Ax = m.atSurfX(j);
        mx[j*(nx+2)] = props.rho * Ax * u[j*(nx+1)];
        mx[j*(nx+2)+nx+1] = props.rho * Ax * u[j*(nx+1)+nx];
        for(int i = 1; i < nx+1; i++) {
            double u1 = u[j*(nx+1)+i-1];
            double u2 = u[j*(nx+1)+i];
            mx[j*(nx+2)+i] = 0.5 * props.rho * (u1 + u2) * Ax;
        }
    }

    my = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    for(int j = 0; j < ny+1; j++) {
        my[j*(nx+1)] = props.rho * m.atSemiSurfY(1,0) * v[j*(nx+2)];
        my[j*(nx+1)+nx] = props.rho * m.atSemiSurfY(nx,1) * v[j*(nx+2)+nx+1];
    }
    for(int i = 1; i < nx; i++) {
        for(int j = 0; j < ny+1; j++) {
            double Ay_left = m.atSemiSurfY(i,1);
            double Ay_right = m.atSemiSurfY(i+1,0);
            double v_left = v[j*(nx+2)+i];
            double v_right = v[j*(nx+2)+i+1];
            my[j*(nx+1)+i] = props.rho * (Ay_left * v_left + Ay_right * v_right);
        }
    }

}

void computeMassFlowsStaggY(double* mx, double* my, const NCMesh m, const double* u, const double* v, const Properties props) {

    int nx = m.getNX();
    int ny = m.getNY();

    mx = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    for(int j = 1; j < ny; j++) {
        for(int i = 0; i < nx+1; i++) {
            double Ax_down = m.atSemiSurfX(j,1);
            double Ax_up = m.atSemiSurfX(j+1,0);
            double u_down = u[j*(nx+1)+i];
            double u_up = u[(j+1)*(nx+1)+i];
            mx[j*(nx+1)+i] = props.rho * (Ax_down * u_down + Ax_up * u_up);
        }
    }


    my = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    for(int i = 1; i < nx+1; i++) {
        double Ay = m.atSurfY(i);
        for(int j = 1; j < ny+1; j++) {
            double v1 = v[(j-1)*(nx+2)+i];
            double v2 = v[j*(nx+2)+i];
            my[j*(nx+2)+i] = 0.5 * props.rho * (v1 + v2) * Ay;
        }
    }
}

void computeVelocitiesStaggX() {

}

void computeVelocitiesStaggY() {
    
}


void computeRu(double* Ru, const NCMesh m, const double* u, const double* v, const Properties props) {

    int nx = m.getNX();
    int ny = m.getNY();

    for(int j = 1; j < ny+1; j++) {
        for(int i = 1; i < nx; i++) {
            // Nodal velocities
            int node = j * (nx + 1) + i;
            double uP = u[node];
            double uW = u[node-1];
            double uE = u[node+1];
            double uS = u[node-(nx+1)];
            double uN = u[node+(nx+1)];
            // Faces velocities
            double uw = schemeCDS(uP, uW);
            double ue = schemeCDS(uP, uE);
            double us = (j > 1 ? 0.5*(uP + uS) : uS);
            double un = (j < ny ? 0.5*(uP + uN) : uN);
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

void computeRv(double* Rv, const NCMesh m, const double* u, const double* v, const Properties props) {

    int nx = m.getNX();
    int ny = m.getNY();

    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny; j++) {
            // Nodal velocities
            int node = j * (nx + 2) + i;
            double vP = v[node];
            double vW = v[node-1];
            double vE = v[node+1];
            double vS = v[node-(nx+2)];
            double vN = v[node+(nx+2)];
            // Face velocities
            double vw = (i > 1 ? 0.5*(vP + vW) : vW);
            double ve = (i < nx ? 0.5*(vP + vE) : vE);
            double vs = schemeCDS(vP, vS);
            double vn = schemeCDS(vP, vN);
            // Areas
            double Ax = m.atSurfX_StaggY(j);
            double Ax_up = m.atSemiSurfX(j+1,0);
            double Ax_down = m.atSemiSurfX(j,1);
            double Ay = m.atSurfY(i);
            // Mass flows
            double mw = props.rho * (u[j*(nx+1)+i-1] * Ax_down + u[(j+1)*(nx+1)+i-1] * Ax_up);
            double me = props.rho * (u[j*(nx+1)+i] * Ax_down + u[(j+1)*(nx+1)+i] * Ax_up);
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

void computePredictorVelocityU(double* u_pred, const int nx, const int ny, const double* u, const double* Ru, const double* Ru_prev, const Properties props, const double tstep) {
    for(int i = 1; i < nx; i++) {
        for(int j = 1; j < ny+1; j++) {
            int node = j * (nx + 1) + i;
            u_pred[node] = u[node] + (1.5 * Ru[node] - 0.5 * Ru_prev[node]) * tstep / props.rho;
        }
    }
}


void computePredictorVelocityV(double* v_pred, const int nx, const int ny, const double* v, const double* Rv, const double* Rv_prev, const Properties props, const double tstep) {
    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny; j++) {
            int node = j * (nx + 2) + i;
            v_pred[node] = v[node] + (1.5 * Rv[node] - 0.5 * Rv_prev[node]) * tstep / props.rho;
        }
    }
}

void computeDiscretizationCoefficients(double* A, double* b, const NCMesh m, const double* u_pred, const double* v_pred, const Properties props, const double tstep) {

    // Mesh sizes
    int nx = m.getNX(); // X-axis control volume count
    int ny = m.getNY(); // Y-axis control volume count

    // Compute discretization coefficients
    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny+1; j++) {
            // Surfaces and distances
            double Ax = m.atSurfX(j);       // X-surface
            double Ay = m.atSurfY(i);       // Y-surface
            double dPE = m.atDistX(i);      // Current node - East node distance
            double dPW = m.atDistX(i-1);    // Current node - West node distance
            double dPN = m.atDistY(j);      // Current node - North node distance
            double dPS = m.atDistY(j-1);    // Current node - South node distance
            // Discretization coefficients
            int node = j * (nx + 2) + i;    // Node number
            A[5*node] = Ay / dPS;           // South node
            A[5*node+1] = Ax / dPW;         // West node
            A[5*node+2] = Ax / dPE;         // East node
            A[5*node+3] = Ay / dPN;         // North node
            A[5*node+4] = A[5*node] + A[5*node+1] + A[5*node+2] + A[5*node+3];  // Current node
            // Independent term
            int node_east = j * (nx + 1) + i;           // X-staggered east node number
            int node_west = j * (nx + 1) + i - 1;       // X-staggered west node number
            int node_north = j * (nx + 2) + i;          // Y-staggered north node number
            int node_south = (j - 1) * (nx + 2) + i;    // Y-staggered south node number
            b[node] = (-props.rho / tstep) * (u_pred[node_east] * Ax - u_pred[node_west] * Ax + v_pred[node_north] * Ay - v_pred[node_south] * Ay);    // Independent term
        }
    }
}

void computeVelocity(double* u, double* v, const NCMesh m, const double* u_pred, const double* v_pred, const double* p, const Properties props, const double tstep, double& maxDiff) {

    // Mesh size
    int nx = m.getNX(); // Number of control volumes along x axis
    int ny = m.getNY(); // Number of control volumes along y axis

    // X-component of velocity
    for(int i = 1; i < nx; i++) {
        for(int j = 1; j < ny+1; j++) {
            int node = j * (nx + 2) + i;                                // Pressure node number
            double px = (p[node+1] - p[node]) / m.atDistX(i);           // Partial derivative of pressure with respect to x
            node = j * (nx + 1) + i;                                    // X-velocity node number
            double u_next = u_pred[node] - (tstep / props.rho) * px;    // X-velocity at time n+1
            maxDiff = std::max(maxDiff, std::abs(u[node] - u_next));    // Update maximum difference
            u[node] = u_next;                                           // Update velocity
        }
    }

    // Y-component of velocity
    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny; j++) {
            int node = j * (nx + 2) + i;                                // Pressure node number
            double py = (p[node+(nx+2)] - p[node]) / m.atDistY(j);      // Partial derivative of pressure with respect to y
            node = j * (nx + 2) + i;                                    // Y-velocity node number
            double v_next = v_pred[node] - (tstep / props.rho) * py;    // Y-velocity at time instant n+1
            maxDiff = std::max(maxDiff, std::abs(v[node] - v_next));    // Update maximum difference to check convergence
            v[node] = v_next;                                           // Update velocity
        }
    }
}

void computeVelocityU(double* u, const NCMesh m, const double* u_pred, const double* p, const Properties props, const double tstep, double& maxDiff) {

    // Mesh size
    int nx = m.getNX(); // Number of control volumes along x axis
    int ny = m.getNY(); // Number of control volumes along y axis

    // X-component of velocity
    for(int i = 1; i < nx; i++) {
        for(int j = 1; j < ny+1; j++) {
            int node = j * (nx + 2) + i;                                // Pressure node number
            double px = (p[node+1] - p[node]) / m.atDistX(i);           // Partial derivative of pressure with respect to x
            node = j * (nx + 1) + i;                                    // X-velocity node number
            double u_next = u_pred[node] - (tstep / props.rho) * px;    // X-velocity at time n+1
            maxDiff = std::max(maxDiff, std::abs(u[node] - u_next));    // Update maximum difference
            u[node] = u_next;                                           // Update velocity
        }
    }

}

void computeVelocityV(double* v, const NCMesh m, const double* v_pred, const double* p, const Properties props, const double tstep, double& maxDiff) {

    // Mesh size
    int nx = m.getNX(); // Number of control volumes along x axis
    int ny = m.getNY(); // Number of control volumes along y axis



    // Y-component of velocity
    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny; j++) {
            int node = j * (nx + 2) + i;                                // Pressure node number
            double py = (p[node+(nx+2)] - p[node]) / m.atDistY(j);      // Partial derivative of pressure with respect to y
            node = j * (nx + 2) + i;                                    // Y-velocity node number
            double v_next = v_pred[node] - (tstep / props.rho) * py;    // Y-velocity at time instant n+1
            maxDiff = std::max(maxDiff, std::abs(v[node] - v_next));    // Update maximum difference to check convergence
            v[node] = v_next;                                           // Update velocity
        }
    }

}

void computeTimeStep(double &tstep, const NCMesh m, const double* u, const double* v, const Properties props, const double tol) {

    // Mesh size
    int nx = m.getNX(); // Number of control volumes along x axis
    int ny = m.getNY(); // Number of control volumes along y axis

    // Compute convective time limit
    double tc = DBL_MAX;    // Large initial value
    // Convective time limit for x axis
    for(int i = 0; i < nx+1; i++) {
        for(int j = 0; j < ny+2; j++) {
            double vel = std::abs(u[j*(nx+1)+i]);
            if(vel > tol)
                tc = std::min(tc, m.atDistX(i) / vel);
        }
    }
    // Convective time limit for y axis
    for(int i = 0; i < nx+2; i++) {
        for(int j = 0; j < ny+1; j++) {
            double vel = std::abs(v[j*(nx+2)+i]);
            if(vel > tol)
                tc = std::min(tc, m.atDistY(j) / vel);
        }
    }

    // Diffusive time limit
    double td = DBL_MAX;    // Large initial value
    double nu = props.mu / props.rho;   // Kinematic viscosity
    // Diffuse time limit for x axis
    for(int i = 0; i < nx+1; i++)
        td = std::min(td, m.atDistX(i) * m.atDistX(i) / nu);
    // Diffuse time limit for y axis
    for(int j = 0; j < ny+1; j++)
        td = std::min(td, m.atDistY(j) * m.atDistY(j) / nu);
    td *= 0.5;

    // Compute new time step
    tstep = 0.2 * std::min(tc, td);
}


void updateOperatorR(double* Ru_prev, double* Rv_prev, const double* Ru, const double* Rv, const int nx, const int ny) {
    // Update Ru
    for(int i = 1; i < nx; i++) {
        for(int j = 1; j < ny+1; j++) {
            int node = j * (nx + 1) + i;
            Ru_prev[node] = Ru[node];
        }
    }
    // Update Rv
    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny; j++) {
            int node = j * (nx + 2) + i;
            Rv_prev[node] = Rv[node];
        }
    }
}

void computeCenteredNodesVelocities(double* u_col, double* v_col, const int nx, const int ny, const double* u, const double* v) {

    // COMPUTE U_COL
    // x-velocity at the center of the non-staggered control volumes
    // Left and right columns
    for(int j = 0; j < ny+2; j++) {
        u_col[j*(nx+2)] = u[j*(nx+1)];
        u_col[j*(nx+2)+nx+1] = u[j*(nx+1)+nx];
    }
    // Central columns
    for(int i = 1; i < nx+1; i++) {
        for(int j = 0; j < ny+2; j++) {
            double uW = u[j*(nx+1)+i-1];
            double uE = u[j*(nx+1)+i];
            u_col[j*(nx+2)+i] = schemeCDS(uW, uE);
        }
    }

    // Compute V_COL
    // y-velocity at the center of the non-staggered control volumes
    // Lower and upper rows
    for(int i = 0; i < nx+2; i++) {
        v_col[i] = v[i];
        v_col[(ny+1)*(nx+2)+i] = v[ny*(nx+2)+i];
    }
    // Central rows
    for(int j = 1; j < ny+1; j++) {
        for(int i = 0; i < nx+2; i++) {
            double vN = v[j*(nx+2)+i];
            double vS = v[(j-1)*(nx+2)+i];
            v_col[j*(nx+2)+i] = schemeCDS(vS, vN);
        }
    }


}

void printVelocityToFile(const NCMesh m, double* u_col, double* v_col, const char* filename, const int precision) {

    std::ofstream file;
    file.open(filename);

    if(!file.is_open()) {
        printf("Error. Could not open file\n");
        file.close();
        return;
    }

    int nx = m.getNX();
    int ny = m.getNY();

    file << std::setprecision(precision) << std::fixed;
    for(int i = 0; i < nx+2; i++) {
        for(int j = 0; j < ny+2; j++) {
            int node = j * (nx + 2) + i;
            double norm = std::sqrt(u_col[node] * u_col[node] + v_col[node] * v_col[node]);
            file << m.atNodeX(i) << " " << m.atNodeY(j) << " " << norm << std::endl;
        }
        file << std::endl;
    }


    file.close();

}

void printVelocityUToFile(const NCMesh m, double* u_col, const char* filename, const int precision) {

    std::ofstream file;
    file.open(filename);

    if(!file.is_open()) {
        printf("Error. Could not open file\n");
        file.close();
        return;
    }

    int nx = m.getNX();
    int ny = m.getNY();

    file << std::setprecision(precision) << std::fixed;
    for(int i = 0; i < nx+2; i++) {
        for(int j = 0; j < ny+2; j++) {
            int node = j * (nx + 2) + i;
            file << m.atNodeX(i) << " " << m.atNodeY(j) << " " << u_col[node] << std::endl;
        }
        file << std::endl;
    }


    file.close();

}

void printVelocityVToFile(const NCMesh m, double* v_col, const char* filename, const int precision) {

    std::ofstream file;
    file.open(filename);

    if(!file.is_open()) {
        printf("Error. Could not open file\n");
        file.close();
        return;
    }

    int nx = m.getNX();
    int ny = m.getNY();

    file << std::setprecision(precision) << std::fixed;
    for(int i = 0; i < nx+2; i++) {
        for(int j = 0; j < ny+2; j++) {
            int node = j * (nx + 2) + i;
            file << m.atNodeX(i) << " " << m.atNodeY(j) << " " << v_col[node] << std::endl;
        }
        file << std::endl;
    }


    file.close();
}


void printPressureToFile(const NCMesh m, double* p, const char* filename, const int precision) {

    std::ofstream file;
    file.open(filename);

    if(!file.is_open()) {
        printf("Error. Could not open file\n");
        file.close();
        return;
    }

    int nx = m.getNX();
    int ny = m.getNY();

    file << std::setprecision(precision) << std::fixed;
    for(int i = 0; i < nx+2; i++) {
        for(int j = 0; j < ny+2; j++) {
            int node = j * (nx + 2) + i;
            file << m.atNodeX(i) << " " << m.atNodeY(j) << " " << p[node] << std::endl;
        }
        file << std::endl;
    }


    file.close();
}

// Lid-driven cavity

void lid_driven::mainLoop(const double rho, const double mu, const double u_ref, const double p_ref, const double L, const int nx, const int ny, const double tstep0, const int maxIt, const double tol, const double sstol) {

    // General variables
    Properties props = {rho, mu};   // Fluid properties
    NCMesh m(L, L, 1, nx, ny);      // Node centered staggered mesh
    double t = 0;           // Current time
    int it = 0;             // Iteration counter
    bool steady = false;    // Steady state bool. true = steady state reached, false = steady state not reached
    double tstep = tstep0;  // Initial time step

    // Fluid variables
    double* u;  // X-velocity array
    double* v;  // Y-velocity array
    double* p;  // Pressure array
    allocateFluidVariables(nx, ny, u, v, p);                // Allocate the previous arrays
    lid_driven::setInitialMaps(u, v, p, m, u_ref, p_ref);   // Set initial values

    // Operators for the fractional step method
    double* Ru;         // Operator Ru at time n
    double* Rv;         // Operator Rv at time n
    double* Ru_prev;    // Operator Ru at time n-1
    double* Rv_prev;    // Operator Rv at time n-1
    allocateOperatorR(nx, ny, Ru, Rv, Ru_prev, Rv_prev);    // Allocate the previous arrays
    computeRu(Ru_prev, m, u, v, props);                     // Compute initial value of Ru
    computeRv(Rv_prev, m, u, v, props);                     // Compute initial value of Rv

    // Predictor velocities
    double* u_pred;     // X-component of predictor velocity
    double* v_pred;     // Y-component of predictor velocity
    allocatePredictorVelocities(nx, ny, u_pred, v_pred);    // Allocate the previous arrays

    // Linear system variables
    double* A;          // Linear system matrix
    double* b;          // Linear system vector
    allocateLinearSystemVariables(nx, ny, A, b);    // Allocate the previous arrays

    std::chrono::steady_clock::time_point begin, end;

    // Run the loop until steady state is reached
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
            printf("Error solving the linear system. Convergence not achieved.\n");
            return;
        }
        // Compute the velocities at time n+1
        double maxDerivative = 0;   // Max derivative of velocity with respect to time
        computeVelocity(u, v, m, u_pred, v_pred, p, props, tstep, maxDerivative);
        maxDerivative /= tstep;
        // Update current time and number of iterations
        t += tstep;
        it++;
        // Time elapsed
        end = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1e6;


        printf("%6d %5s %10.5f %5s %10.5f %5s %10.5e %5s %10d %5s %.3f\n", it, "", t, "", tstep, "", maxDerivative, "", exitCodeGS, "", elapsed);

        if(it % 10 == 0) {

            double* u_col = (double*) calloc((nx+2)*(ny+2), sizeof(double));
            double* v_col = (double*) calloc((nx+2)*(ny+2), sizeof(double));
            if(!u_col) {
                printf("Error: could not allocate enough memory for u_col\n");
                return;
            }
            if(!v_col) {
                printf("Error: could not allocate enough memory for v_col\n");
                return;
            }
            computeCenteredNodesVelocities(u_col, v_col, nx, ny, u, v);

            int Re = std::floor(props.rho * u_ref * L / props.mu);
            std::string filename = "../plots/vel_" + std::to_string(nx) + "_" + std::to_string(ny) + "_" + std::to_string(Re) + ".txt";
            printVelocityToFile(m, u_col, v_col, filename.c_str(), 5);
            filename = "../plots/u" + std::to_string(nx) + "_" + std::to_string(ny) + "_" + std::to_string(Re) + ".txt";
            printVelocityUToFile(m, u_col, "sim/u.txt", 5);
            filename = "../plots/v" + std::to_string(nx) + "_" + std::to_string(ny) + "_" + std::to_string(Re) + ".txt";
            printVelocityVToFile(m, v_col, "sim/v.txt", 5);
            filename = "../plots/p" + std::to_string(nx) + "_" + std::to_string(ny) + "_" + std::to_string(Re) + ".txt";
            printPressureToFile(m, p, "sim/p.txt", 5);

            free(u_col);
            free(v_col);

        }

        // Check steady state condition
        if(maxDerivative < sstol)
            steady = true;
        else {
            // Compute next time step
            computeTimeStep(tstep, m, u, v, props, tol);
            // Update operator R
            std::memcpy(Ru_prev, Ru, (nx+1)*(ny+2)*sizeof(double));
            std::memcpy(Rv_prev, Rv, (nx+2)*(ny+1)*sizeof(double));
        }
    }

    // Free arrays
    free(u);
    free(v);
    free(p);
    free(Ru);
    free(Rv);
    free(Ru_prev);
    free(Rv_prev);
    free(u_pred);
    free(v_pred);
    free(A);
    free(b);
}

void lid_driven::setInitialMaps(double* u, double* v, double* p, const NCMesh m, const double u_ref, const double p_ref) {

    // Mesh size
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


void lid_driven::setBoundaryPredictorVelocities(double* u_pred, double* v_pred, const int nx, const int ny, const double u_ref) {

    // PREDICTOR VELOCITY X
    // Lower and upper boundaries
    for(int i = 0; i < nx+1; i++) {
        u_pred[i] = 0;              // Lower boundary
        int j = ny + 1;
        u_pred[j*(nx+1)+i] = u_ref; // Upper boundary
    }

    // Left and right boundaries
    for(int j = 1; j < ny+1; j++) {
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

    // Mesh size
    int nx = m.getNX(); // Number of control volumes along x axis
    int ny = m.getNY(); // Number of control volumes along y axis

    // lOWER AND UPPER BOUNDARIES (including corner nodes, i.e. singular nodes)
    for(int i = 0; i < nx+2; i++) {
        // Lower boundary
        int node = i;       // Node number
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
