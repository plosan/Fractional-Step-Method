#include "NCMesh.h"

NCMesh::NCMesh() : built(false), lx(0), ly(0), lz(0), nx(0), ny(0) {

}

NCMesh::NCMesh(double _lx, double _ly, double _lz, int _nx, int _ny) : built(false), lx(_lx), ly(_ly), lz(_lz), nx(_nx), ny(_ny) {

    // Tells whether or not a step of the construction of the mesh could be carried out successfully
    //  - true: the step was finished successfully
    //  - false: the step could not be finished since calloc could not allocate enough memory for some member variable
    bool construction = computeFaceXY();
    if(!construction) {
        printf("Mesh: error. The mesh will halt its construction\n");
        return;
    }

    construction = computeNodeXY();
    if(!construction) {
        printf("Mesh: error. The mesh will halt its construction\n");
        return;
    }

    construction = computeDistXY();
    if(!construction) {
        printf("Mesh: error. The mesh will halt its construction\n");
        return;
    }

    construction = computeSurfXY();
    if(!construction) {
        printf("Mesh: error. The mesh will halt its construction\n");
        return;
    }

    construction = computeVol();
    if(!construction) {
        printf("Mesh: error. The mesh will halt its construction\n");
        return;
    }

    built = true;
    printf("Mesh: the mesh was successfully constructed\n");

}

// For a uniform mesh, the function computes:
//  - faceX: X-position of the faces perpendicular to the X axis
//  - faceY: Y-position of the faces perpendicular to the Y axis
// Return value:
//  - false: if the function was unable to allocate memory either for faceX or faceY
//  - true: otherwise
bool NCMesh::computeFaceXY() {
    // Compute faceX: X-position of the faces perpendicular to the X axis in a uniform mesh
    faceX = (double*) calloc(nx+1, sizeof(double));
    if(!faceX) {
        printf("Mesh: error. Could not allocate enough memory for faceX\n");
        return false;
    }
    // faceX[0] = 0;       // First face at x = 0
    // faceX[nx] = lx;     // Last face at x = lx
    // double d = lx/nx;   // X-size of a control volume in a uniform mesh
    // for(int i = 1; i < nx; i++)
    //     faceX[i] = i*d;

    double d = lx/nx;
    for(int i = 0; i < nx+1; i++)
        faceX[i] = i*d;

    // Compute faceY: Y-position of the faces perpendicular to the Y axis in a uniform mesh
    faceY = (double*) calloc(ny+1, sizeof(double));
    if(!faceY) {
        printf("Mesh: error. Could not allocate enough memory for faceY\n");
        return false;
    }
    // faceY[0] = 0;       // First face at y = 0
    // faceY[ny] = ly;     // Last face at y = ly
    // d = ly/ny;          // Y-size of a control volume in a uniform mesh
    // for(int j = 1; j < ny; j++)
    //     faceY[j] = j*d;
    d = ly/ny;
    for(int j = 0; j < ny+1; j++)
        faceY[j] = j*d;

    // Everything good so far
    printf("Mesh: faces position computed successfully\n");
    return true;
}

// For a uniform or non-uniform mesh, the function computes:
//  - nodeX: position of nodes in X coordinate
//  - nodeY: position of nodes in Y coordinate
// Return value:
//  - false: if the function was unable to allocate memory either for nodeX or nodeY
//  - true: otherwise
bool NCMesh::computeNodeXY() {
    // Compute nodeX: X-position of nodes
    nodeX = (double*) calloc(nx+2, sizeof(double));
    if(!nodeX) {
        printf("Mesh: error. Could not allocate enough memory for nodeX\n");
        return false;
    }
    nodeX[0] = 0;       // First node at x = 0
    nodeX[nx+1] = lx;   // Last node at x = lx
    for(int i = 1; i < nx+1; i++)
        nodeX[i] = 0.5*(faceX[i-1] + faceX[i]); // Other nodes equidistant from adjacent faces

    // Compute nodeY: Y-position of nodes
    nodeY = (double*) calloc(ny+2, sizeof(double));
    if(!nodeY) {
        printf("Mesh: error. Could not allocate enough memory for nodeY\n");
        return false;
    }
    nodeY[0] = 0;       // First node at y = 0
    nodeY[ny+1] = ly;   // Last node at y = ly
    for(int j = 1; j < ny+1; j++)
        nodeY[j] = 0.5*(faceY[j-1] + faceY[j]); // Other nodes equisdistant from adjacent faces

    // Everything good so far
    printf("Mesh: nodes position computed successfully\n");
    return true;
}

// For a uniform or non-uniform mesh, the function computes:
//  - distX: distances between nodes in the X coordinate
//  - distY: distances between nodes in the Y coordinate
// Return value:
//  - false: if the function was unable to allocate memory either for distX or distY
//  - true: otherwise
bool NCMesh::computeDistXY() {
    // Compute distX: distances between nodes in the X coordinate
    distX = (double*) calloc(nx+1, sizeof(double));
    if(!distX) {
        printf("Mesh: error. Could not allocate enough memory for distX\n");
        return false;
    }
    for(int i = 0; i < nx+1; i++)
        distX[i] = nodeX[i+1] - nodeX[i];

    // Compute distY: distances between nodes in the Y coordinate
    distY = (double*) calloc(ny+1, sizeof(double));
    if(!distY)  {
        printf("Mesh: error. Could not allocate enough memory for distY\n");
        return false;
    }
    for(int j = 0; j < ny+1; j++)
        distY[j] = nodeY[j+1] - nodeY[j];

    // Everything good so far
    printf("Mesh: distances between nodes computed successfully\n");
    return true;
}

// For a uniform or non-uniform mesh, the function computes:
//  - surfX: surface of the faces perpendicular to the X axis
//  - surfY: surface of the faces perpendicular to the Y axis
// Return value:
//  - false: if the function was unable to allocate memory either for surfX or surfY
//  - true: otherwise
bool NCMesh::computeSurfXY() {
    // Compute surfX: surface of the faces perpendicular to the X axis
    surfX = (double*) calloc(ny+2, sizeof(double));
    if(!surfX) {
        printf("Mesh: error. Could not allocate enough memory for surfX\n");
        return false;
    }
    surfX[0] = 0;
    surfX[ny+1] = 0;
    for(int j = 1; j < ny+1; j++)
        surfX[j] = (faceY[j] - faceY[j-1]) * lz;

    // Compute surfY: surface of the faces perpendicular to the Y axis
    surfY = (double*) calloc(nx+2, sizeof(double));
    if(!surfY) {
        printf("Mesh: error. Could not allocate enough memory for surfY\n");
        return false;
    }
    surfY[0] = 0;
    surfY[nx+1] = 0;
    for(int i = 1; i < nx+1; i++)
        surfY[i] = (faceX[i] - faceX[i-1]) * lz;

    // Everything good so far
    printf("Mesh: surfaces computed successfully\n");
    return true;
}

// For a uniform or non-uniform mesh, the function computes:
//  - vol: volume of the control volumes
// Return value:
//  - false: if the function was unable to allocate memory for vol
//  - true: otherwise
bool NCMesh::computeVol() {
    // Compute vol: volume of the control volumes
    vol = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    if(!vol) {
        printf("Mesh: error. Could not allocate enough memory for vol\n");
        return false;
    }
    for(int i = 0; i < nx+2; i++)
        for(int j = 0; j < ny+2; j++)
            vol[i+j*(nx+2)] = surfX[j] * surfY[i] / lz;

    // Everything good so far
    printf("Mesh: volumes computed successfully\n");
    return true;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GETTERS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Returns built
bool NCMesh::isBuilt() const {
    return built;
}

// Returns lx
double NCMesh::getLX() const {
    if(!built) {
        printf("Mesh: error. The mesh is not built. Returning 0 for lx\n");
        return 0;
    }
    return lx;
}

// Returns ly
double NCMesh::getLY() const {
    if(!built) {
        printf("Mesh: error. The mesh is not built. Returning 0 for lY\n");
        return 0;
    }
    return ly;
}

// Returns lz
double NCMesh::getLZ() const {
    if(!built) {
        printf("Mesh: error. The mesh is not built. Returning 0 for lZ\n");
        return 0;
    }
    return lz;
}

// Returns nx
int NCMesh::getNX() const {
    if(!built) {
        printf("Mesh: error. The mesh is not built. Returning 0 for nx\n");
        return 0;
    }
    return nx;
}

// Returns ny
int NCMesh::getNY() const {
    if(!built) {
        printf("Mesh: error. The mesh is not built. Returning 0 for ny\n");
        return 0;
    }
    return ny;
}

// Returns nodeX
double* NCMesh::getNodeX() const {
    if(!built) {
        printf("Mesh: error. The mesh is not built. Returning nullptr for nodeX\n");
        return nullptr;
    }
    return nodeX;
}

// Returns nodeY
double* NCMesh::getNodeY() const {
    if(!built) {
        printf("Mesh: error. The mesh is not built. Returning nullptr for nodeY\n");
        return nullptr;
    }
    return nodeY;
}

// Returns distX
double* NCMesh::getDistX() const {
    if(!built) {
        printf("Mesh: error. The mesh is not built. Returning nullptr for distX\n");
        return nullptr;
    }
    return distX;
}

// Returns distY
double* NCMesh::getDistY() const {
    if(!built) {
        printf("Mesh: error. The mesh is not built. Returning nullptr for distY\n");
        return nullptr;
    }
    return distY;
}

// Returns faceX
double* NCMesh::getFaceX() const {
    if(!built) {
        printf("Mesh: error. The mesh is not built. Returning nullptr for faceX\n");
        return nullptr;
    }
    return faceX;
}

// Returns faceY
double* NCMesh::getFaceY() const {
    if(!built) {
        printf("Mesh: error. The mesh is not built. Returning nullptr for faceY\n");
        return nullptr;
    }
    return faceY;
}

// Returns surfX
double* NCMesh::getSurfX() const {
    if(!built) {
        printf("Mesh: error. The mesh is not built. Returning nullptr for surfX\n");
        return nullptr;
    }
    return surfX;
}

// Returns surfY
double* NCMesh::getSurfY() const {
    if(!built) {
        printf("Mesh: error. The mesh is not built. Returning nullptr for surfY\n");
        return nullptr;
    }
    return surfY;
}

// Returns vol
double* NCMesh::getVol() const {
    if(!built) {
        printf("Mesh: error. The mesh is not built. Returning nullptr for vol\n");
        return nullptr;
    }
    return vol;
}

// Returns nodeX[i] unsafely (does not check if 0 <= i < nx+2)
double NCMesh::atNodeX(int i) const {
    return nodeX[i];
}

// Returns nodeY[j] unsafely (does not check if 0 <= j < ny+2)
double NCMesh::atNodeY(int j) const {
    return nodeY[j];
}

// Returns distX[i] unsafely (does not check if 0 <= i < nx+1)
double NCMesh::atDistX(int i) const {
    return distX[i];
}

// Returns distY[j] unsafely (does not check if 0 <= j < ny+1)
double NCMesh::atDistY(int j) const {
    return distY[j];
}

// Returns faceX[i] unsafely (does not check if 0 <= i < nx+1)
double NCMesh::atFaceX(int i) const {
    return faceX[i];
}

// Returns faceY[j] unsafely (does not check if 0 <= j < ny+1)
double NCMesh::atFaceY(int j) const {
    return faceY[j];
}

// Returns surfX[j] unsafely (does not check if 0 <= j < ny+2)
double NCMesh::atSurfX(int j) const {
    return surfX[j];
}

// Returns surfY[i] unsafely (does not check if 0 <= i < nx+2)
double NCMesh::atSurfY(int i) const {
    return surfY[i];
}

// Returns vol[j*nx+i] unsafely (does not check if 0 <= i < nx+2 and 0 <= j < ny+2)
double NCMesh::atVol(int i, int j) const {
    return vol[i+j*(nx+2)];
}


// Displays the following mesh member variables:
//  - built
//  - lx, ly, lz
//  - nx, ny
//  - Node number and location
//  - distX, distY
//  - surfX, surfY
//  - vol
void NCMesh::printMeshData() const {

    // GENERAL DATA 1
    printf("\nGeneral data\n");
    printf("%10s%2s|%2s%s\n\n", "Built", "", "", (built ? "Yes" : "No"));

    if(built) {

        // GENERAL DATA 2
        printf("%10s%3s%10s%10s%10s\n", "", "", "x", "y", "z");
        printf("%10s%2s|%10.5f%10.5f%10.5f\n", "Length", "", lx, ly, lz);
        printf("%10s%2s|%10d%10d%10d\n", "CVs", "", nx, ny, 1);

        // NODE NUMBER AND LOCATION
        printf("\nNode number and location\n");
        // Print column numbers
        printf("%13s", "");
        for(int i = 0; i < nx+2; i++)
            printf("%5d", i);
        // Print horizontal line separating column numbers and table content
        printf("\n%12s", "");
        for(int i = 0; i < 5*(nx+2)+1; i++)
            printf("-");
        printf("\n");
        // Print node number and location
        for(int j = ny+1; j >= 0; j--) {
            printf("%10d%2s|", j, "");
            for(int i = 0; i < nx+2; i++)
                printf("%5d", i+j*(nx+2));
            printf("\n");
        }

        // Print distX
        printf("\nDistance X\n");

        printf("%10s%2s|", "i", "");
        for(int i = 0; i < nx+1; i++)
            printf("%10d", i);

        printf("\n%10s%2s|", "distX[i]", "");
        for(int i = 0; i < nx+1; i++)
            printf("%10.5f", distX[i]);
        printf("\n");

        // Print distY
        printf("\nDistance Y\n");
        printf("%10s%5s%s\n", "j", "", "distY[j]");
        printf("%7s", "");
        for(int i = 0; i < 16; i++)
            printf("-");
        printf("\n");
        for(int j = ny; j >= 0; j--)
            printf("%10d%5s%.5f\n", j, "", distY[j]);

        // SURFACE X
        printf("\nSurfaces X\n");
        // Print table header
        printf("%10s%5s%s\n", "j", "", "surfX[j]");
        printf("%7s", "");
        for(int i = 0; i < 16; i++)
            printf("-");
        printf("\n");
        // Print surfaces X
        for(int j = ny+1; j >= 0; j--)
            printf("%10d%5s%.5f\n", j, "", surfX[j]);

        // SURFACE Y
        printf("\nSurfaces Y\n");
        // Print column numbers
        printf("%10s%2s|", "i", "");
        for(int i = 0; i < nx+2; i++)
            printf("%10d", i);
        // Print surfaces Y
        printf("\n%10s%2s|", "surfY[i]", "");
        for(int i = 0; i < nx+2; i++)
            printf("%10.5f", surfY[i]);
        printf("\n");

        // VOLUMES ASSOCIATED TO EVERY NODE
        printf("\nVolumes * 1e2\n");
        // Print column numbers
        printf("%8s", "");
        for(int i = 0; i < nx+2; i++)
            printf("%10d", i);
        printf("\n%7s", "");
        // Print horizontal line separating column numbers and table content
        for(int i = 0; i < 10*(nx+2)+1; i++)
            printf("-");
        printf("\n");

        for(int j = ny+1; j >= 0; j--) {
            printf("%5d%2s%s", j, "", "|");
            for(int i = 0; i < nx+2; i++)
                printf("%10.2f", 1e2*vol[i+j*(nx+2)]);
            printf("\n");
        }
    }
}

// Saves to files some mesh parameters that can be plotted later on using plotmesh.gnu. Parameters and files:
//  - lx, ly        File: data/domain_data.dat
//  - nodeX, nodeY  File: data/node_data.dat
//  - faceX         File: data/facex_data.dat
//  - faceY         File: data/facey_data.dat
void NCMesh::saveMeshData() const {

    if(!built) {
        printf("Mesh: Error. The mesh is not built. Unable to save mesh data\n");
        return;
    }

    FILE *fp;
    // Save domain size
    fp = fopen("data/domain_size.dat", "w");
    if(!fp) {
        printf("Mesh: error. Could not open file 'domain_data.dat'\n");
        return;
    }
    fprintf(fp, "%.3f %.3f %.3f %.3f\n", 0.0, 0.0, lx, 0.0);
    fprintf(fp, "%.3f %.3f %.3f %.3f\n", 0.0, 0.0, 0.0, ly);
    fprintf(fp, "%.3f %.3f %.3f %.3f\n", lx, 0.0, lx, ly);
    fprintf(fp, "%.3f %.3f %.3f %.3f\n", 0.0, ly, lx, ly);
    fclose(fp);

    // Save nodeX and nodeY
    fp = fopen("data/node_data.dat", "w");
    if(!fp) {
        printf("Mesh: error. Could not open file 'node_data.dat'\n");
        return;
    }
    for(int i = 0; i < nx+2; i++)
        for(int j = 0; j < ny+2; j++)
            fprintf(fp, "%.3f %.3f\n", nodeX[i], nodeY[j]);
    fclose(fp);

    // Save faceX location
    fp = fopen("data/facex_data.dat", "w");
    if(!fp) {
        printf("Mesh: error. Could not open file 'facex_data.dat'\n");
        return;
    }
    for(int i = 0; i < nx+1; i++)
        fprintf(fp, "%.3f %.3f %.3f %.3f\n", faceX[i], 0.0, faceX[i], ly);
    fclose(fp);

    // Save faceY location
    fp = fopen("data/facey_data.dat", "w");
    if(!fp) {
        printf("Mesh: error. Could not open file 'facey_data.dat'\n");
        return;
    }
    for(int j = 0; j < ny+1; j++)
        fprintf(fp, "%.3f %.3f %.3f %.3f\n", 0.0, faceY[j], lx, faceY[j]);
    fclose(fp);

}
