
#include "RCGrid.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RCGrid::RCGrid() : built(false), lx(0), ly(0), lz(0), nx(0), ny(0) {

}

RCGrid::RCGrid(double _lx, double _ly, double _lz, int _nx, int _ny) :
    built(false), lx(_lx), ly(_ly), lz(_lz), nx(_nx), ny(_ny) {

    bool error = computeNodeXY();
    if(!error) {
        printf("Error: the mesh will stop its construction\n");
        return;
    }

    error = computeDistXY();
    if(!error) {
        printf("Error: the mesh will stop its construction\n");
        return;
    }

    error = computeFaceXY();
    if(!error) {
        printf("Error: the mesh will stop its construction\n");
        return;
    }

    error = computeSurfXY();
    if(!error) {
        printf("Error: the mesh will stop its construction\n");
        return;
    }

    error = computeVol();
    if(!error) {
        printf("Error: the mesh will stop its construction\n");
        return;
    }

    built = true;
    printf("The mesh was successfully constructed\n");
}

RCGrid::RCGrid(double _lx, double _ly, double _lz, int _nx, int _ny, double* _nodeX, double* _nodeY) {

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS TO COMPUTE GRID PARAMETERS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Compute the position of nodes in X and Y coordinates for a uniform mesh
bool RCGrid::computeNodeXY() {
    // Compute nodeX, position of nodes in X coordinate for a uniform mesh
    nodeX = (double*) calloc(nx, sizeof(double));
    if(!nodeX) {
        printf("Error: could not allocate enough memory for nodeX\n");
        return false;
    }
    nodeX[0] = 0;
    nodeX[nx-1] = lx;
    double d = lx/(nx - 1);
    for(int i = 1; i < nx-1; i++)
        nodeX[i] = i*d;

    // Compute nodeY, position of nodes in Y coordinate for a uniform mesh
    nodeY = (double*) calloc(ny, sizeof(double));
    if(!nodeY) {
        printf("Error: could not allocate enough memory for nodeY\n");
        return false;
    }
    nodeY[0] = 0;
    nodeY[ny-1] = ly;
    d = ly/(ny - 1);
    for(int j = 1; j < ny-1; j++)
        nodeY[j] = j*d;

    // Everything good so far
    return true;
}

// Compute the distances between nodes in X coordinate and Y coordinate
bool RCGrid::computeDistXY() {
    // Compute distX, distance between nodes in X coordinate
    distX = (double*) calloc(nx-1, sizeof(double));
    if(!distX) {
        printf("Error: could not allocate enough memory for distX\n");
        return false;
    }
    for(int i = 0; i < nx-1; i++)
        distX[i] = nodeX[i+1] - nodeX[i];

    // Compute distY, distance between nodes in Y coordinate
    distY = (double*) calloc(ny-1, sizeof(double));
    if(!distY) {
        printf("Error: could not allocate enough memory for distY\n");
        return false;
    }
    for(int j = 0; j < ny-1; j++)
        distY[j] = nodeY[j+1] - nodeY[j];

    // Everything good so far
    return true;
}

// Compute the location of faces perpendicular to the X axis and Y axis
bool RCGrid::computeFaceXY() {
    // Compute faceX, location of the faces perpendicular to the X axis
    faceX = (double*) calloc(nx+1, sizeof(double));
    if(!faceX) {
        printf("Error: could not allocate enough memory for faceX\n");
        return false;
    }
    faceX[0] = 0;
    faceX[nx] = lx;
    for(int i = 1; i < nx; i++)
        faceX[i] = 0.5*(nodeX[i-1] + nodeX[i]);

    // Compute faceY, location of the faces perpendicular to the Y axis
    faceY = (double*) calloc(ny+1, sizeof(double));
    if(!faceY) {
        printf("Error: could not allocate enough memory for faceY\n");
        return false;
    }
    faceY[0] = 0;
    faceY[ny] = ly;
    for(int j = 1; j < ny; j++)
        faceY[j] = 0.5*(nodeY[j-1] + nodeY[j]);

    // Everything good so far
    return true;
}

// Compute the surface of the faces perpendicular to the X axis and Y axis
bool RCGrid::computeSurfXY() {
    // Compute surfX, surface of the faces perpendicular to the X axis
    surfX = (double*) calloc(ny, sizeof(double));
    if(!surfX) {
        printf("Error: could not allocate enough memory for surfX\n");
        return false;
    }
    for(int j = 0; j < ny; j++)
        surfX[j] = (faceY[j+1] - faceY[j])*lz;

    // Compute surfY, surface of the faces perpendicular to the Y axis
    surfY = (double*) calloc(nx, sizeof(double));
    if(!surfY) {
        printf("Error: could not allocate enough memory for surfY\n");
        return false;
    }
    for(int i = 0; i < nx; i++)
        surfY[i] = (faceX[i+1] - faceX[i])*lz;

    // Everything good so far
    return true;
}

// Compute the volume of each control volume
bool RCGrid::computeVol() {
    // Compute vol, volume of each control volume
    vol = (double*) calloc(nx*ny, sizeof(double));
    if(!vol) {
        printf("Error: could not allocate enough memory for vol\n");
        return false;
    }
    for(int j = 0; j < ny; j++)
        for(int i = 0; i < nx; i++)
            vol[j*nx+i] = surfX[j] * surfY[i] * lz;

    // Everything good so far
    return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GETTERS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Returns built
bool RCGrid::isBuilt() const {
    return built;
}

// Returns lx
double RCGrid::getLX() const {
    if(!built) {
        printf("Error: the mesh is not built. Returning 0 for lx\n");
        return 0;
    }
    return lx;
}

// Returns ly
double RCGrid::getLY() const {
    if(!built) {
        printf("Error: the mesh is not built. Returning 0 for ly\n");
        return 0;
    }
    return ly;
}

// Returns lz
double RCGrid::getLZ() const {
    if(!built) {
        printf("Error: the mesh is not built. Returning 0 for lz\n");
        return 0;
    }
    return lz;
}

// Returns nx
int RCGrid::getNX() const {
    if(!built) {
        printf("Error: the mesh is not built. Returning 0 for nx\n");
        return 0;
    }
    return nx;
}

// Returns nu
int RCGrid::getNY() const {
    if(!built) {
        printf("Error: the mesh is not built. Returning 0 for ny\n");
        return 0;
    }
    return ny;
}

// Returns nodeX
double* RCGrid::getNodeX() const {
    if(!built) {
        printf("Error: the mesh is not built. Returning nullptr for nodeX\n");
        return nullptr;
    }
    return nodeX;
}

// Returns nodeY
double* RCGrid::getNodeY() const {
    if(!built) {
        printf("Error: the mesh is not built. Returning nullptr for nodeY\n");
        return nullptr;
    }
    return nodeY;
}

// Returns distX
double* RCGrid::getDistX() const {
    if(!built) {
        printf("Error: the mesh is not built. Returning nullptr for distX\n");
        return nullptr;
    }
    return distX;
}

// Returns distY
double* RCGrid::getDistY() const {
    if(!built) {
        printf("Error: the mesh is not built. Returning nullptr for distY\n");
        return nullptr;
    }
   return distY;
}

// Returns faceX
double* RCGrid::getFaceX() const {
    if(!built) {
        printf("Error: the mesh is not built. Returning nullptr for faceX\n");
        return nullptr;
    }
    return faceX;
}

// Returns faceY
double* RCGrid::getFaceY() const {
    if(!built) {
        printf("Error: the mesh is not built. Returning nullptr for faceY\n");
        return nullptr;
    }
    return faceY;
}

// Returns surfX
double* RCGrid::getSurfX() const {
    if(!built) {
        printf("Error: the mesh is not built. Returning nullptr for surfX\n");
        return nullptr;
    }
    return surfX;
}

// Returns surfY
double* RCGrid::getSurfY() const {
    if(!built) {
        printf("Error: the mesh is not built. Returning nullptr for surfY\n");
        return nullptr;
    }
    return surfY;
}

// Returns vol
double* RCGrid::getVol() const {
    if(!built) {
        printf("Error: the mesh it not built. Returning nullptr for vol\n");
        return nullptr;
    }
    return vol;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SAFE GETTERS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Returns nodeX[i] safely (checks if the object is built and if 0 <= i < nx)
double RCGrid::getNodeX(int i) const {
    if(!built) {
        printf("Error: the mesh is not built. Returning 0 for nodeX[%d]\n", i);
        return 0;
    }
    if(i < 0 || i >= nx) {
        printf("Error: index i = %d out of bounds (0 - %d). Returning 0 for nodeX[%d]\n", i, nx-1, i);
        return 0;
    }
    return nodeX[i];
}

// Returns nodeY[j] safely (checks if the object is built and if 0 <= j < ny)
double RCGrid::getNodeY(int j) const {
    if(!built) {
        printf("Error: the mesh is not built. Returning 0 for nodeY[%d]\n", j);
        return 0;
    }
    if(j < 0 || j >= ny) {
        printf("Error: index j = %d out of bounds (0 - %d). Returning 0 for nodeY[%d]\n", j, ny-1, j);
        return 0;
    }
    return nodeY[j];
}

// Returns distX[i] safely (checks if the object is built and if 0 <= i < nx-1)
double RCGrid::getDistX(int i) const {
    if(!built) {
        printf("Error: the mesh is not built. Returning 0 for distX[%d]\n", i);
        return 0;
    }
    if(i < 0 || i >= nx-1) {
        printf("Error: index i = %d out of bounds (0 - %d). Returning 0 for distX[%d]\n", i, nx-2, i);
        return 0;
    }
    return distX[i];
}

// Returns distY[j] safely (checks if the object is built and if 0 <= j < ny-1)
double RCGrid::getDistY(int j) const {
    if(!built) {
        printf("Error: the mesh is not built. Returning 0 for distY[%d]\n", j);
        return 0;
    }
    if(j < 0 || j >= ny-1) {
        printf("Error: index j = %d out of bounds (0 - %d). Returning 0 for distY[%d]\n", j, ny-2, j);
        return 0;
    }
    return distY[j];
}

// Returns faceX[i] safely (checks if the object is built and if 0 <= i < nx+1)
double RCGrid::getFaceX(int i) const {
    if(!built) {
        printf("Error: the mesh is not built. Returning 0 for faceX[%d]\n", i);
        return 0;
    }
    if(i < 0 || i >= nx+1) {
        printf("Error: index i = %d out of bounds (0 - %d). Returning 0 for faceX[%d]\n", i, nx, i);
        return 0;
    }
    return faceX[i];
}

// Returns faceY[j] safely (checks if the object is built and if 0 <= j < ny+1)
double RCGrid::getFaceY(int j) const {
    if(!built) {
        printf("Error: the mesh is not built. Returning 0 for faceY[%d]\n", j);
        return 0;
    }
    if(j < 0 || j >= ny+1) {
        printf("Error: index j = %d out of bounds (0 - %d). Returning 0 for faceY[%d]\n", j, ny, j);
        return 0;
    }
    return faceY[j];
}

// Returns surfX[j] safely (checks if the object is built and if 0 <= j < ny)
double RCGrid::getSurfX(int j) const {
    if(!built) {
        printf("Error: the mesh is not built. Returning 0 for surfX[%d]\n", j);
        return 0;
    }
    if(j < 0 || j >= ny) {
        printf("Error: index j = %d out of bounds (0 - %d). Returning 0 for surfX[%d]\n", j, ny-1, j);
        return 0;
    }
    return surfX[j];
}

// Returns surfY[i] safely (checks if the object is built and if 0 <= i < nx)
double RCGrid::getSurfY(int i) const {
    if(!built) {
        printf("Error: the mesh is not built. Returning 0 for surfY[%d]\n", i);
        return 0;
    }
    if(i < 0 || i >= nx) {
        printf("Error: index i = %d out of bounds (0 - %d). Returning 0 for surfY[%d]\n", i, nx-1, i);
        return 0;
    }
    return surfY[i];
}

// Returns vol[i+j*nx] safely (checks if the object is built and if 0 <= i < nx nad 0 <= j < ny)
double RCGrid::getVol(int i, int j) const {
    if(!built) {
        printf("Error: the mesh is not built. Returning 0 for vol[%d,%d]\n", i, j);
        return 0;
    }
    if(i < 0 || i >= nx) {
        printf("Error: index i = %d out of bounds (0 - %d). Returning 0 for vol[%d][%d]\n", i, nx-1, i, j);
        return 0;
    }
    if(j < 0 || j >= ny) {
        printf("Error: index j = %d out of bounds (0 - %d). Returning 0 for vol[%d][%d]\n", j, ny-1, i, j);
        return 0;
    }
    return vol[i+j*nx];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UNSAFE GETTERS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Returns nodeX[i] unsafely
double RCGrid::atNodeX(int i) const {
    return nodeX[i];
}

// Returns nodeY[j] unsafely
double RCGrid::atNodeY(int j) const {
    return nodeY[j];
}

// Returns distX[i] unsafely
double RCGrid::atDistX(int i) const {
    return distX[i];
}

// Returns distY[j] unsafely
double RCGrid::atDistY(int j) const {
    return distY[j];
}

// Returns faceX[i] unsafely
double RCGrid::atFaceX(int i) const {
    return faceX[i];
}

// Returns faceY[j] unsafely
double RCGrid::atFaceY(int j) const {
    return faceY[j];
}

// Returns surfX[j] unsafely
double RCGrid::atSurfX(int j) const {
    return surfX[j];
}

// Returns surfY[i] unsafely
double RCGrid::atSurfY(int i) const {
    return surfY[i];
}

// Returns vol[i+j*nx] unsafely
double RCGrid::atVol(int i, int j) const {
    return vol[i+j*nx];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// INFORMATION FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Displays the following mesh member variables:
//  - built
//  - lx, ly, lz
//  - nx, ny
//  - Node number and location
//  - distX, distY
//  - surfX, surfY
//  - vol
void RCGrid::printMeshData() const {

    // General data
    printf("\nGeneral data\n");
    printf("%10s%2s|%2s%s\n\n", "Built", "", "", (built ? "Yes" : "No"));

    if(built) {

        printf("%10s%3s%10s%10s%10s\n", "", "", "x", "y", "z");
        printf("%10s%2s|%10.5f%10.5f%10.5f\n", "Length", "", lx, ly, lz);
        printf("%10s%2s|%10d%10d%10d\n", "Nodes", "", nx, ny, 1);

        // Print node id and location
        printf("\nNode number and location\n");
        printf("%13s", "");
        for(int i = 0; i < nx; i++)
            printf("%5d", i);
        printf("\n%12s", "");
        for(int i = 0; i < 5*nx+1; i++)
            printf("-");
        printf("\n");

        for(int j = 0; j < ny; j++) {
            printf("%10d%2s|", ny-1-j, "");
            for(int i = 0; i < nx; i++)
                printf("%5d", (ny-1-j)*nx+i);
            printf("\n");
        }

        // Print distX
        printf("\nDistance X\n");

        printf("%10s%2s|", "i", "");
        for(int i = 0; i < nx-1; i++)
            printf("%10d", i);

        printf("\n%10s%2s|", "distX[i]", "");
        for(int i = 0; i < nx-1; i++)
            printf("%10.5f", distX[i]);
        printf("\n");

        // Print distY
        printf("\nDistance Y\n");
        printf("%10s%5s%s\n", "j", "", "distY[j]");
        for(int j = ny-2; j >= 0; j--)
            printf("%10d%5s%.5f\n", j, "", distY[j]);

        // Surface X
        printf("\nSurfaces X\n");
        printf("%10s%5s%s\n", "j", "", "surfX[j]");
        for(int j = ny-1; j >= 0; j--)
            printf("%10d%5s%.5f\n", j, "", surfX[j]);

        // Surface Y
        printf("\nSurfaces Y\n");
        printf("%10s%2s|", "i", "");
        for(int i = 0; i < nx; i++)
            printf("%10d", i);
        printf("\n");
        printf("%10s%2s|", "surfY[i]", "");
        for(int i = 0; i < nx; i++)
            printf("%10.5f", surfY[i]);
        printf("\n");

        // Volumes
        printf("\nVolumes\n");
        printf("%8s", "");
        for(int i = 0; i < nx; i++)
            printf("%10d", i);
        printf("\n%7s", "");
        for(int i = 0; i < 10*nx+1; i++)
            printf("-");
        printf("\n");
        for(int j = ny-1; j >= 0; j--) {
            printf("%5d%2s%s", j, "", "|");
            for(int i = 0; i < nx; i++)
                printf("%10.5f", vol[j*nx+i]);
            printf("\n");
        }
    }
}

// Saves to files some mesh parameters that can be plotted later on using plotmesh.gnu. Parameters and files:
//  - lx, ly        File: data/domain_data.dat
//  - nodeX, nodeY  File: data/node_data.dat
//  - faceX         File: data/facex_data.dat
//  - faceY         File: data/facey_data.dat
void RCGrid::saveMeshData() const {

    if(!built) {
        printf("Error: the mesh is not built. Unable to save mesh data\n");
        return;
    }

    FILE *fp;   // Pointer to file
    // Save domain size
    fp = fopen("data/domain_data.dat", "w");
    if(!fp) {
        printf("Error: could not open file 'domain_data.dat'\n");
        return;
    }
    fprintf(fp, "%.3f %.3f %.3f %.3f\n", 0.0, 0.0, lx, 0.0);
    fprintf(fp, "%.3f %.3f %.3f %.3f\n", 0.0, 0.0, 0.0, ly);
    fprintf(fp, "%.3f %.3f %.3f %.3f\n", lx, 0.0, lx, ly);
    fprintf(fp, "%.3f %.3f %.3f %.3f\n", 0.0, ly, lx, ly);
    fclose(fp);

    // Save node location
    fp = fopen("data/node_data.dat", "w");
    if(!fp) {
        printf("Error: could not open file 'node_data.dat'\n");
        return;
    }
    for(int i = 0; i < nx; i++)
        for(int j = 0; j < ny; j++)
            fprintf(fp, "%.3f %.3f\n", nodeX[i], nodeY[j]);
    fclose(fp);

    // Save faceX location
    fp = fopen("data/facex_data.dat", "w");
    if(!fp) {
        printf("Error: could not open file 'facex_data.dat'\n");
        return;
    }
    for(int i = 0; i < nx+1; i++)
        fprintf(fp, "%.3f %.3f %.3f %.3f\n", faceX[i], 0.0, faceX[i], ly);
    fclose(fp);

    // Save faceY location
    fp = fopen("data/facey_data.dat", "w");
    if(!fp) {
        printf("Error: could not open file 'facey_data.dat'\n");
        return;
    }
    for(int j = 0; j < ny+1; j++)
        fprintf(fp, "%.3f %.3f %.3f %.3f\n", 0.0, faceY[j], lx, faceY[j]);
    fclose(fp);
}
