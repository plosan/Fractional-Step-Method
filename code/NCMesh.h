#ifndef NCMESH_H
#define NCMESH_H

#include <iostream>

class NCMesh {

private:
    bool built;     // Whether the mesh is built (functions) or not. Used to tell if getters should access member variables
    double lx;      // Length of the domain in the X coordinate
    double ly;      // Length of the domain in the Y coordinate
    double lz;      // Length of the domain in the Z coordinate
    int nx;         // Number of control volumes for discretisation in the X coordinate
    int ny;         // Number of control volumes for discretisation in the Y coordinate
    double* nodeX;  // Position of nodes in the X coordinate. Size: nx+2
    double* nodeY;  // Position of nodes in the Y coordinate. Size: ny+2
    double* distX;  // Distances between nodes in the X coordinate. Size: nx+1
    double* distY;  // Distances between nodes in the Y coordinate. Size: ny+1
    double* faceX;  // Position in the X coordinate of the faces perpendicular to the X axis. Size: nx+1
    double* faceY;  // Position in the Y coordinate of the faces perpendicular to the Y axis. Size: ny+1
    double* surfX;  // Surface of the faces perpendicular to the X axis. One surface associated to each node along the X axis. Size: ny+2
    double* surfY;  // Surface of the faces perpendicular to the Y axis. One surface associated to each node along the Y axis. Size: nx+2
    double* vol;    // Volume of the control volumes. For simplicity, one control volume associated to each node although some have 0 volume. Size: (nx+2)*(ny+2)

    bool computeFaceXY();
    bool computeNodeXY();
    bool computeDistXY();
    bool computeSurfXY();
    bool computeVol();

    bool computeSurfXY2();
    bool computeVol2();

public:

    // Constructors
    NCMesh();
    NCMesh(double _lx, double _ly, double _lz, int _nx, int _ny);
    NCMesh(double _lx, double _ly, double _lz, int _nx, int _ny, double* _faceX, double* _faceY);

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
    double atNodeX(int) const;     // Returns nodeX[i] unsafely (does not check if 0 <= i < nx+2)
    double atNodeY(int) const;     // Returns nodeY[j] unsafely (does not check if 0 <= j < ny+2)
    double atDistX(int) const;     // Returns distX[i] unsafely (does not check if 0 <= i < nx+1)
    double atDistY(int) const;     // Returns distY[j] unsafely (does not check if 0 <= j < ny+1)
    double atFaceX(int) const;     // Returns faceX[i] unsafely (does not check if 0 <= i < nx+1)
    double atFaceY(int) const;     // Returns faceY[j] unsafely (does not check if 0 <= j < ny+1)
    double atSurfX(int) const;     // Returns surfX[j] unsafely (does not check if 0 <= j < ny+2)
    double atSurfY(int) const;     // Returns surfY[i] unsafely (does not check if 0 <= i < nx+2)
    double atVol(int, int) const;  // Returns vol[j*nx+i] unsafely (does not check if 0 <= i < nx+2 and 0 <= j < ny+2)

    // Info functions
    void printMeshData() const; // Prints the mesh parameters
    void saveMeshData() const;  // Saves the mesh parameters to different files to plot it later on


};

#endif /* NCMESH_H */
