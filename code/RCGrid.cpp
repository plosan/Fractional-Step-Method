
#include "RCGrid.h"

RCGrid::RCGrid(double _lx, double _ly, double _lz, int _nx, int _ny) :
    lx(_lx), ly(_ly), lz(_lz), nx(_nx), ny(_ny) {
    computeNodeXY();
    computeFaceXY();
    computeSurfXY();
    computeVol();

    for(int j = ny-1; j >= 0; j--) {
        for(int i = 0; i < nx; i++)
            printf("%5d", j*nx+i);
        printf("\n");
    }
}

RCGrid::RCGrid(double _lx, double _ly, double _lz, int _nx, int _ny, double* _nodeX, double* _nodeY) {

}

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

    // Compute nodeOY, position of nodes in Y coordinate
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

void RCGrid::printMeshData() const {

    printf("%15s%10s%10s%10s\n", "", "x", "y", "z");
    printf("%15s%10.3f%10.3f%10.3f\n", "Length [m]", lx, ly, lz);
    printf("%15s%10d%10d%10d\n", "Nodes", nx, ny, 1);

    // Print node id and location
    printf("\nNode number and location\n");
    printf("%8s", "");
    for(int i = 0; i < nx; i++)
        printf("%5d", i);
    printf("\n%8s", "");
    for(int i = 0; i < 5*nx; i++)
        printf("-");
    printf("\n");

    for(int j = 0; j < ny; j++) {
        printf("%5d%2s%s", ny-1-j, "", "|");
        for(int i = 0; i < nx; i++)
            printf("%5d", (ny-1-j)*nx+i);
        printf("\n");
    }

    // Surface X
    printf("\nSurfaces X\n");
    printf("%5s%5s%s\n", "j", "", "surfX[j] [m2]");
    for(int j = 0; j < ny; j++)
        printf("%5d%5s%.3f\n", j, "", surfX[j]);

    // Surface Y
    printf("\nSurfaces Y\n");
    printf("%5s%5s%s\n", "i", "", "surfY[i] [m2]");
    for(int i = 0; i < nx; i++)
        printf("%5d%5s%.3f\n", i, "", surfY[i]);

    // Volumes
    printf("\nVolumes\n");
    for(int j = ny-1; j >= 0; j--) {
        printf("%5d%2s%s", j, "", "|");
        for(int i = 0; i < nx; i++)
            printf("%10.5f", vol[j*nx+i]);
        printf("\n");
    }
}

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
