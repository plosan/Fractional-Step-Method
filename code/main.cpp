#include <iostream>
#include <vector>
#include <cstring>
#include "RCGrid.h"


int main(int argc, char* argv[]) {

    RCGrid m(1, 1, 1, 7, 5);

    m.saveMeshData();
    m.printMeshData();


}
