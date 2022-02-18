
#include "RCGrid.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RCGrid::RCGrid(double _lx, double _ly, double _lz, int _nx, int _ny) :
    lx(_lx), ly(_ly), lz(_lz), nx(_nx), ny(_ny) {
    computeNodeXY();
    computeDistXY();
    computeFaceXY();
    computeSurfXY();
    computeVol();
}

RCGrid::RCGrid(double _lx, double _ly, double _lz, int _nx, int _ny, double* _nodeX, double* _nodeY) {

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS TO COMPUTE GRID PARAMETERS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Compute the position of nodes in X and Y coordinates
void RCGrid::computeNodeXY() {
    // Compute nodeX, position of nodes in X coordinate
    nodeX = (double*) calloc(nx, sizeof(double));
    if(!nodeX) {
        printf("Error: could not allocate enough memory for nodeX\n");
        return;
    }
    nodeX[0] = 0;
    nodeX[nx-1] = lx;
    double d = lx/(nx - 1);
    for(int i = 1; i < nx-1; i++)
        nodeX[i] = i*d;

    // Compute nodeY, position of nodes in Y coordinate
    nodeY = (double*) calloc(ny, sizeof(double));
    if(!nodeY) {
        printf("Error: could not allocate enough memory for nodeY\n");
        return;
    }
    nodeY[0] = 0;
    nodeY[ny-1] = ly;
    d = ly/(ny - 1);
    for(int j = 1; j < ny-1; j++)
        nodeY[j] = j*d;
}

// Compute the distances between nodes in X coordinate and Y coordinate
void RCGrid::computeDistXY() {
    // Compute distX, distance between nodes in X coordinate
    distX = (double*) calloc(nx-1, sizeof(double));
    if(!distX) {
        printf("Error: could not allocate enough memory for distX\n");
        return;
    }
    for(int i = 0; i < nx-1; i++)
        distX[i] = nodeX[i+1] - nodeX[i];

    // Compute distY, distance between nodes in Y coordinate
    distY = (double*) calloc(ny-1, sizeof(double));
    if(!distY) {
        printf("Error: could not allocate enough memory for distY\n");
        return;
    }
    for(int j = 0; j < ny-1; j++)
        distY[j] = nodeY[j+1] - nodeY[j];
}

// Compute the location of faces perpendicular to the X axis and Y axis
void RCGrid::computeFaceXY() {
    // Compute faceX, location of the faces perpendicular to the X axis
    faceX = (double*) calloc(nx+1, sizeof(double));
    if(!faceX) {
        printf("Error: could not allocate enough memory for faceX\n");
        return;
    }
    faceX[0] = 0;
    faceX[nx] = lx;
    for(int i = 1; i < nx; i++)
        faceX[i] = 0.5*(nodeX[i-1] + nodeX[i]);

    // Compute faceY, location of the faces perpendicular to the Y axis
    faceY = (double*) calloc(ny+1, sizeof(double));
    if(!faceY) {
        printf("Error: could not allocate enough memory for faceY\n");
        return;
    }
    faceY[0] = 0;
    faceY[ny] = ly;
    for(int j = 1; j < ny; j++)
        faceY[j] = 0.5*(nodeY[j-1] + nodeY[j]);
}

// Compute the surface of the faces perpendicular to the X axis and Y axis
void RCGrid::computeSurfXY() {
    // Compute surfX, surface of the faces perpendicular to the X axis
    surfX = (double*) calloc(ny, sizeof(double));
    if(!surfX) {
        printf("Error: could not allocate enough memory for surfX\n");
        return;
    }
    for(int j = 0; j < ny; j++)
        surfX[j] = (faceY[j+1] - faceY[j])*lz;

    // Compute surfY, surface of the faces perpendicular to the Y axis
    surfY = (double*) calloc(nx, sizeof(double));
    if(!surfY) {
        printf("Error: could not allocate enough memory for surfY\n");
        return;
    }
    for(int i = 0; i < nx; i++)
        surfY[i] = (faceX[i+1] - faceX[i])*lz;
}

// Compute the volume of each control volume
void RCGrid::computeVol() {
    // Compute vol, volume of each control volume
    vol = (double*) calloc(nx*ny, sizeof(double));
    if(!vol) {
        printf("Error: could not allocate enough memory for vol\n");
        return;
    }
    for(int j = 0; j < ny; j++)
        for(int i = 0; i < nx; i++) {
            vol[j*nx+i] = surfX[j]*surfY[i];
        }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GETTERS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int RCGrid::getNX() const {
    return nx;
}

int RCGrid::getNY() const {
    return ny;
}

double* RCGrid::getNodeX() const {
    return nodeX;
}

double* RCGrid::getNodeY() const {
    return nodeY;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// INFORMATION FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Prints the mesh parameters
void RCGrid::printMeshData() const {

    // General data
    printf("\nGeneral data\n");
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

// Saves the mesh parameters to different files to plot it later on
void RCGrid::saveMeshData() const {

    FILE *fp;
    // Save domain size
    fp = fopen("domain_data.dat", "w");
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
    fp = fopen("node_data.dat", "w");
    if(!fp) {
        printf("Error: could not open file 'node_data.dat'\n");
        return;
    }
    for(int i = 0; i < nx; i++)
        for(int j = 0; j < ny; j++)
            fprintf(fp, "%.3f %.3f\n", nodeX[i], nodeY[j]);
    fclose(fp);

    // Save faceX location
    fp = fopen("facex_data.dat", "w");
    if(!fp) {
        printf("Error: could not open file 'facex_data-dat'\n");
        return;
    }
    for(int i = 0; i < nx+1; i++)
        fprintf(fp, "%.3f %.3f %.3f %.3f\n", faceX[i], 0.0, faceX[i], ly);
    fclose(fp);

    // Save faceY location
    fp = fopen("facey_data.dat", "w");
    if(!fp) {
        printf("Error: could not open file 'facey_data.dat'\n");
        return;
    }
    for(int j = 0; j < ny+1; j++)
        fprintf(fp, "%.3f %.3f %.3f %.3f\n", 0.0, faceY[j], lx, faceY[j]);
    fclose(fp);
}
