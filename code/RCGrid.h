#ifndef RCGRID_H
#define RCGRID_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>

class RCGrid {

private:
    bool built;     // Tells whether or not the object is built
    double lx;      // Length of the domain in X coordinate
    double ly;      // Length of the domain in Y coordinate
    double lz;      // Length of the domain in Z coordinate
    int nx;         // Number of nodes in X coordinate
    int ny;         // Number of nods in Y coordinate
    double* nodeX;  // Nodes location in X coordinate. Size: nx
    double* nodeY;  // Nodes location in Y coordinate. Size: ny
    double* distX;  // Distance between nodes in X coordinate. Size: nx-1
    double* distY;  // Distance between nodes in Y coordinate. Size: ny-1
    double* faceX;  // Location of the faces perpendicular to the X axis. Size: nx+1
    double* faceY;  // Location of the faces perpendicular to the Y axis. Size: ny+1
    double* surfX;  // Surface of the faces perpendicular to the X axis. Size: ny
    double* surfY;  // Surface of the faces perpendicular to the Y axis. Size: nx
    double* vol;    // Volume of each control volume. Size: nx*ny

    // Compute grid parameters
    bool computeNodeXY();   // Compute the position of nodes in X and Y coordinates for a uniform mesh
    bool computeDistXY();   // Compute the distances between nodes in X coordinate and Y coordinate
    bool computeFaceXY();   // Compute the location of faces perpendicular to the X axis and Y axis
    bool computeSurfXY();   // Compute the surface of the faces perpendicular to the X axis and Y axis
    bool computeVol();      // Compute the volume of each control volume



public:
    // Constructors
    RCGrid();
    RCGrid(double _lx, double _ly, double _lz, int _nx, int _ny);
    RCGrid(double _lx, double _ly, double _lz, int _nx, int _ny, double* _nodeX, double* _nodeY);

    // Getters
    bool isBuilt() const;       // Returns built
    double getLX() const;       // Returns lx
    double getLY() const;       // Returns ly
    double getLZ() const;       // Returns lz
    int getNX() const;          // Returns nx
    int getNY() const;          // Returns ny
    double* getNodeX() const;   // Returns nodeX
    double* getNodeY() const;   // Returns nodeY
    double* getDistX() const;   // Returns distX
    double* getDistY() const;   // Returns distY
    double* getFaceX() const;   // Returns faceX
    double* getFaceY() const;   // Returns faceY
    double* getSurfX() const;   // Returns surfX
    double* getSurfY() const;   // Returns surfY
    double* getVol() const;     // Returns vol

    // Safe getters
    double getNodeX(int) const;     // Returns nodeX[i] safely (checks if the object is built and if 0 <= i < nx)
    double getNodeY(int) const;     // Returns nodeY[j] safely (checks if the object is built and if 0 <= j < ny)
    double getDistX(int) const;     // Returns distX[i] safely (checks if the object is built and if 0 <= i < nx-1)
    double getDistY(int) const;     // Returns distY[j] safely (checks if the object is built and if 0 <= j < ny-1)
    double getFaceX(int) const;     // Returns faceX[i] safely (checks if the object is built and if 0 <= i < nx+1)
    double getFaceY(int) const;     // Returns faceY[j] safely (checks if the object is built and if 0 <= j < ny+1)
    double getSurfX(int) const;     // Returns surfX[j] safely (checks if the object is built and if 0 <= j < ny)
    double getSurfY(int) const;     // Returns surfY[i] safely (checks if the object is built and if 0 <= i < nx)
    double getVol(int, int) const;  // Returns vol[j*nx+i] safely (checks if the object is built and if 0 <= i < nx and 0 <= j < ny)

    // Unsafe getters
    double atNodeX(int) const;     // Returns nodeX[i] unsafely (does not check if 0 <= i < nx)
    double atNodeY(int) const;     // Returns nodeY[j] unsafely (does not check if 0 <= j < ny)
    double atDistX(int) const;     // Returns distX[i] unsafely (does not check if 0 <= i < nx-1)
    double atDistY(int) const;     // Returns distY[j] unsafely (does not check if 0 <= j < ny-1)
    double atFaceX(int) const;     // Returns faceX[i] unsafely (does not check if 0 <= i < nx+1)
    double atFaceY(int) const;     // Returns faceY[j] unsafely (does not check if 0 <= j < ny+1)
    double atSurfX(int) const;     // Returns surfX[j] unsafely (does not check if 0 <= j < ny)
    double atSurfY(int) const;     // Returns surfY[i] unsafely (does not check if 0 <= i < nx)
    double atVol(int, int) const;  // Returns vol[j*nx+i] unsafely (does not check if 0 <= i < nx and 0 <= j < ny)

    // Info functions
    void printMeshData() const; // Prints the mesh parameters
    void saveMeshData() const;  // Saves the mesh parameters to different files to plot it later on

};

#endif // RCGRID_H
