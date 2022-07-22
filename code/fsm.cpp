#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cfloat>

#include "fsm.h"
#include "schemes.h"
#include "solver.h"


void allocateFluidVariables(const int nx, const int ny, double* &u, double* &v, double* &p) {
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

void computeVelocityCollocatedMesh(double* u_col, double* v_col, const int nx, const int ny, const double* u, const double* v) {

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

// Lid-driven cavity

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
