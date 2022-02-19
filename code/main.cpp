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

    // m.saveMeshData();
    m.printMeshData();


}
