#include <iostream>
#include <vector>
#include <cstring>
#include "RCGrid.h"


int main(int argc, char* argv[]) {

    double lx = 1;
    double ly = 1;
    double lz = 1;
    int nx = 7;
    int ny = 5;
    // RCGrid m;
    RCGrid m(lx, ly, lz, nx, ny);




    printf("%10s = %s\n", "is built", (m.isBuilt() ? "Yes" : "No"));
    printf("%10s = %.5f\n", "lx", m.getLX());
    printf("%10s = %.5f\n", "ly", m.getLY());
    printf("%10s = %.5f\n", "lz", m.getLZ());
    printf("%10s = %d\n", "nx", m.getNX());
    printf("%10s = %d\n", "ny", m.getNY());

    double* f = m.getVol();
    printf("%10s = %p\n", "vol", m.getVol());
    printf("%10s = %ld\n", "size", sizeof(f));
    for(int j = 0; j < m.getNY(); j++) {
        for(int i = 0; i < m.getNX(); i++)
            printf("%5d%5d%2s%.5f\n", i, j, "", f[i+j*m.getNX()]);
        printf("\n");
    }


    // printf("%10s = %p\n", "nodeY", m.getNodeY());
    // printf("%10s = %p\n", "distX", m.getDistX());
    // printf("%10s = %p\n", "distY", m.getDistY());
    // printf("%10s = %p\n", "faceX", m.getFaceX());
    // printf("%10s = %p\n", "faceY", m.getFaceY());
    // printf("%10s = %p\n", "surfX", m.getSurfX());
    // printf("%10s = %p\n", "surfY", m.getSurfY());
    // printf("%10s = %p\n", "vol", m.getVol());
    // m.getNodeX();
    // m.getNodeY();
    // m.getDistX();
    // m.getDistY();
    // m.getFaceX();
    // m.getFaceY();
    // m.getSurfX();
    // m.getSurfY();
    // m.getVol();

    // m.saveMeshData();
    m.printMeshData();


}
