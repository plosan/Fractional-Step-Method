#include <iostream>
#include <vector>
#include <cstring>
#include "RCGrid.h"
#include "NCMesh.h"
#include "schemes.h"

#define NX 10
#define NY 10
#define INDEX(I,J) I+J*(NY+2)

void computeRu();
void computeRv();
void computePredictorVelocityX();
void computePredictorVelocityY();
double schemeCDS(double A, double B);

int main(int argc, char* argv[]) {

    double L = 1;
    int nx = NX;
    int ny = NY;

    NCMesh m(L, L, 1, nx, ny);
    m.saveMeshData();
    m.printMeshData();


    // X-component of velocity
    double* u = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    if(!u) {
        printf("Error: could not allocate enough memory for u\n");
        return false;
    }

    // Y-component of velocity
    double* v = (double*) calloc((nx+2)*(ny+2), sizeof(double));
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


}


double schemeCDS(double A, double B) {
    return 0.5*(A + B);
}

void computeRu() {

    // Horizontal component of operator R
    // For corner nodes ((0,0) and (nx,0) and (0,ny+1) and (nx,ny+1)) Ru is zero
    double* Ru = (double*) calloc((nx+1)*(ny+2), sizeof(double));

    // Left wall
    for(int j = 1; j < ny+1; j++) {
        // Velocity ue
        double uP = u[j*(nx+1)];
        double uE = u[1+j*(nx+1)];
        double ue = schemeCDS(uP, uE);
        // Velocity un
        double un = 0;
        if(j < ny) {
            un = schemeCDS(u[j*(nx+1)], u[(j+1)*(nx+1)]);
        }
        // Velocity us
        double us = 0;
        if(j > 1) {
            us = schemeCDS(u[j*(nx+1)], u[(j-1)*(nx+1)]);
        }
        // Velocities us_A and us_B
        double usA = 0;
        double usB = v[1+(j-1)*(nx+2)];
        double AsB = distX[0]*lz;
        // Velocities un_A, un_B and areas
        double unA = 0;
        double unB = v[1+j*(nx+2)];
        double AnB = distX[0]*lz;
        // Mass flow mn
        double mn = rho*unB*AnB;
        // Mass flow me
        double me = 0.5*rho*(u[0+j*(nx+1)] + u[1+j*(nx+1)])*surfX[j];
        // Mass flow ms
        double ms = rho*usB*AsB;
        // Mass flow mw
        double mw = 0;

        // Velocity un
        double uP = u[j*(nx+1)];
        double uE = u[1+j*(nx+1)];
        double uS = u[(j-1)*(nx+1)];
        double uN = u[(j+1)*(nx+1)];

        double Aw = distX[0] * lz;
        double An = Aw;
        Ru[j*(nx+1)] = -(me*ue - mw*uw + mn*un - ms*us);
        Ru[j*(nx+1)] += mu*(uE-uP)*surfX[j]/distX[0] + mu*(uN-uP)*An/distY[j] - mu*(uP-uS)*As/distY[j-1];
        Ru[j*(nx+1)] /= volStaggX[j*(nx+1)];
    }

}


// void computePredictorVelocityX(double* up, const RCGrid m, const double* un)
