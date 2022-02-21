#include "NCMesh.h"

NCMesh::NCMesh(double _lx, double _ly, double _lz, int _nx, int _ny) : lx(_lx), ly(_ly), lz(_lz), nx(_nx), ny(_ny) {


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
        printf("Error: could not allocate enough memory for faceX\n");
        return false;
    }
    faceX[0] = 0;       // First face at x = 0
    faceX[nx] = lx;     // Last face at x = lx
    double d = lx/nx;   // X-size of a control volume in a uniform mesh
    for(int i = 1; i < nx; i++)
        faceX[i] = i*d;

    // Compute faceY: Y-position of the faces perpendicular to the Y axis in a uniform mesh
    faceY = (double*) calloc(ny+1, sizeof(double));
    if(!faceY) {
        printf("Error: could not allocate enough memory for faceY\n");
        return false;
    }
    faceY[0] = 0;       // First face at y = 0
    faceY[ny] = ly;     // Last face at y = ly
    d = ly/ny;          // Y-size of a control volume in a uniform mesh
    for(int j = 1; j < ny; j++)
        faceY[j] = j*d;

    // Everything good so far
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
        printf("Error: could not allocate enough memory for nodeX\n");
        return false;
    }
    nodeX[0] = 0;       // First node at x = 0
    nodeX[nx+1] = lx;   // Last node at x = lx
    for(int i = 1; i < nx+1; i++)
        nodeX[i] = 0.5*(faceX[i-1] + faceX[i]); // Other nodes equidistant from adjacent faces

    // Compute nodeY: Y-position of nodes
    nodeY = (double*) calloc(ny+2, sizeof(double));
    if(!nodeY) {
        printf("Error: could not allocate enough memory for nodeY\n");
        return false;
    }
    nodeY[0] = 0;       // First node at y = 0
    nodeY[ny+1] = ly;   // Last node at y = ly
    for(int j = 1; j < ny+1; j++)
        nodeY[j] = 0.5*(faceY[j-1] + faceY[j]); // Other nodes equisdistant from adjacent faces

    // Everything good so far
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
    surfX = (double*) calloc(ny, sizeof(double));
    if(!surfX) {
        printf("Error: could not allocate enough memory for surfX\n");
        return false;
    }
    for(int j = 0; j < ny; j++)
        surfX[j] = (faceY[j+1] - faceY[j]) * lz;  // Surface of the j-th face perpendicular to the X axis

    // Compute surfY: surface of the faces perpendicular to the Y axis
    surfY = (double*) calloc(nx, sizeof(double));
    if(!surfY) {
        printf("Error: could not allocate enough memory for surfY\n");
        return false;
    }
    for(int i = 0; i < nx; i++)
        surfY[i] = (faceX[i+1] - faceX[i]) * lz;

    // Everything good so far
    return true;
}

// For a uniform or non-uniform mesh, the function computes:
//  - vol: volume of the control volumes
// Return value:
//  - false: if the function was unable to allocate memory for vol
//  - true: otherwise
bool NCMesh::computeVol() {
    // Compute vol: volume of the control volumes
    vol = (double*) calloc(nx*ny, sizeof(double));
    if(!vol) {
        printf("Error: could not allocate enough memory for vol\n");
        return false;
    }
    for(int i = 0; i < nx; i++)
        for(int j = 0; j < ny; j++)
            vol[i+j*nx] = surfX[j] * surfY[i] * lz;

    // Everything good so far
    return true;
}
