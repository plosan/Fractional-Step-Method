#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>

class RCGrid {

private:
    double lx;      // Length of the domain in X coordinate
    double ly;      // Length of the domain in Y coordinate
    double lz;      // Length of the domain in Z coordinate
    int nx;         // Number of nodes in X coordinate
    int ny;         // Number of nods in Y coordinate
    double* nodeX;  // Nodes location in X coordinate. Size: nx
    double* nodeY;  // Nodes location in Y coordinate. Size: ny
    double* faceX;  // Location of the faces perpendicular to the X axis. Size: nx+1
    double* faceY;  // Location of the faces perpendicular to the Y axis. Size: ny+1
    double* surfX;  // Surface of the faces perpendicular to the X axis. Size: ny
    double* surfY;  // Surface of the faces perpendicular to the Y axis. Size: nx
    double* vol;    // Volume of each control volume. Size: nx*ny

    void computeNodeXY();
    void computeFaceXY();
    void computeSurfXY();
    void computeVol();


public:
    // Constructors
    RCGrid();
    RCGrid(double _lx, double _ly, double _lz, int _nx, int _ny);
    RCGrid(double _lx, double _ly, double _lz, int _nx, int _ny, double* _nodeX, double* _nodeY);

    // Info functions
    void printMeshData() const;
    void saveMeshData() const;

};
