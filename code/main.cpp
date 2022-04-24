#include <iostream>
#include <vector>
#include <cstring>
#include "RCGrid.h"
#include "NCMesh.h"
#include "schemes.h"

#define NX 10
#define NY 10
#define INDEX(I,J) I+J*(NY+2)

void computePredictorVelocityX();
void computePredictorVelocityY();

int main(int argc, char* argv[]) {

    double L = 1;
    int nx = NX;
    int ny = NY;

    NCMesh m(L, L, 1, nx, ny);
    m.saveMeshData();
    m.printMeshData();

    // RCGrid m(L, L, 1, nx, ny);

    // double* u = (double*) calloc(nx*ny, sizeof(double*));
    // double* v = (double*) calloc(nx*ny, sizeof(double*));

    // m.saveMeshData();
    // m.printMeshData();

    // for(int j = NY+1; j >= 0; j--) {
    //     for(int i = 0; i < NX+2; i++)
    //         printf("(%2d,%2d) %3d %5s", i, j, INDEX(i,j), "");
    //     printf("\n");
    // }


}

bool f(double lx, double ly, double lz, int nx, int ny) {

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



    return true;

}

// void computePredictorVelocityX(double* up, const RCGrid m, const double* un)
