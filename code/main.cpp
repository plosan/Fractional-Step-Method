#include <iostream>
#include <vector>
#include <cstring>
#include "RCGrid.h"
#include "NCMesh.h"
#include "schemes.h"

#define NX 10
#define NY 10
#define INDEX(I,J) I+J*(NY+2)

double* allocateDoubleArray(int n, int m);
void allocateMatrices(double* &u, double* &v, double* &Ru, double* &Rv, const NCMesh m);

void computeRu();
void computeRv();
void computePredictorVelocityX();
void computePredictorVelocityY();

double schemeUDS(double A, double B, double mf);
double schemeCDS(double A, double B);

struct Properties {
    double rho;
    double mu;
};

int main(int argc, char* argv[]) {

    double rho = 0.9982;    // Water density at 20 ºC           [kg/m^3]
    double mu = 1.0016e-3;  // Water dynamic viscosity at 20 ºC [Pa s]
    double uref = 1;        // X-velocity boundary condition    [m/s]
    double pref = 1e5;      // Pressure

    double L = 1;

    int nx = NX;
    int ny = NY;

    NCMesh m(L, L, 1, nx, ny);
    m.saveMeshData();
    m.printMeshData();


    Properties props = {rho, mu};


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

    // Initial maps
    for(int i = 1; i < nx; i++) {
        int j = ny + 1;
        int id = j * (nx + 1) + i;
        u[id] = uref;
    }

    // Operator Ru
    double* Ru = (double*) calloc((nx+1)*(ny+2), sizeof(double));
    if(!u) {
        printf("Error: could not allocate enough memory for Ru\n");
        return false;
    }

    // Try to compute Ru
    double* mx_Ru = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    if(!mx_Ru) {
        printf("Error: could not allocate enough memory for mx_Ru\n");
        return false;
    }
    for(int j = 1; j < ny+1; j++) {

        double Ax = m.getSurfX(j);  // Surface perpendicular to the X axis
        double mw = 0.5 * props.rho * (u[j*(nx+1)] + u[j*(nx+1)+1]) * Ax;
        double me = 0.5*props.rho*(u[j*(nx+1)+1] + u[j*(nx+1)+2]) * Ax;
        for(int i = 1; i < nx; i++) {

            // double uw = schemeUDS(u[j*(nx+1)+i-1], u[j*(nx+1)+i], mw);
            // double ue = schemeUDS(u[j*(nx+1)+i], u[j*(nx+1)+i+1], me);
            double uw = schemeCDS(u[j*(nx+1)+i-1], u[j*(nx+1)+i]);
            double ue = schemeCDS(u[j*(nx+1)+i], u[j*(nx+1)+i+1]);
            double us = schemeCDS(u[(j-1)*(nx+1)+i], u[j*(nx+1)+i]);    // If and only if j > 1
            double un = schemeCDS(u[j*(nx+1)+i], u[(j+1)*(nx+1)+i]);    // If and only if j < ny

            double mn = props.rho*(v[j*(nx+1)+i] * m.getSemiSurfY(i,0) + v[j*(nx+1)+i+1] * m.getSemiSurfY(i,1));
            double ms = props.rho*(v[(j-1)*(nx+1)+i] * m.getSemiSurfY(i,0) + v[(j-1)*(nx+1)+i+1] * m.getSemiSurfY(i,1));

            double uP = u[j*(nx+1)+i];
            double uW = u[j*(nx+1)+i-1];
            double uE = u[j*(nx+1)+i+1];
            double uS = u[(j-1)*(nx+1)+i];
            double uN = u[(j+1)*(nx+1)+i];

            double integral1 = -(me * ue - mw * uw + mn * un - ms * us);
            double integral2 = m.getSurfX(j) * (uE - uP) / m.getDistX(i) - m.getSurfX(j) * (uP - uW) / m.getDistX(i-1);
            integral2 += m.getSurfY(i) * (uN - uP) / m.getDistY(j) - m.getSurfY(i) * (uP - uS) / m.getDistY(j-1);
            integral2 *= props.mu;

            Ru[j*(nx+1)+i] = (integral1 + integral2) / m.atVolStaggX(i,j);


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

// void computeRu() {
//
//     // Horizontal component of operator R
//     // For corner nodes ((0,0) and (nx,0) and (0,ny+1) and (nx,ny+1)) Ru is zero
//     double* u = (double*) calloc((nx+1)*(ny+2), sizeof(double));
//     double* Ru = (double*) calloc((nx+1)*(ny+2), sizeof(double));
//
//
//
// }


// void computePredictorVelocityX(double* up, const RCGrid m, const double* un)
