#include <iostream>
#include <vector>
#include <cstring>
#include "RCGrid.h"
#include "NCMesh.h"
#include "schemes.h"


void computePredictorVelocityX();
void computePredictorVelocityY();

int main(int argc, char* argv[]) {

    double L = 1;
    int nx = 10;
    int ny = 10;

    NCMesh m(L, L, 1, nx, ny);
    m.saveMeshData();
    m.printMeshData();

    // RCGrid m(L, L, 1, nx, ny);

    // double* u = (double*) calloc(nx*ny, sizeof(double*));
    // double* v = (double*) calloc(nx*ny, sizeof(double*));

    // m.saveMeshData();
    // m.printMeshData();


}

// void computePredictorVelocityX(double* up, const RCGrid m, const double* un)
