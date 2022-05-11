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

    double* distFaceX;  // Distances between faces in the X coordinate. Size: nx
    double* distFaceY;  // Distances between faces in the Y coordinate. Size: ny

    double* volStaggX;  // Volume associated to X-staggered nodes. Size: (nx+1)*(ny+2)
    double* volStaggY;  // Volume associated to Y-staggered nodes. Size: (nx+2)*(ny+1)

    bool computeFaceXY();   // Computes faceX and faceY for uniform and non-uniform meshes
    bool computeNodeXY();   // Computes nodeX and nodeY
    bool computeDistXY();   // Computes distX and distY
    bool computeSurfXY();   // Computes surfX and surfY
    bool computeVol();      // Computes vol

    bool computeDistFaceXY();   // Computes distFaceX and distFaceY
    bool computeVolStaggX();    // Computes volStaggX
    bool computeVolStaggY();    // Computes volStaggY

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
    double* getDistFaceX() const;   // Returns distFaceX
    double* getDistFaceY() const;   // Returns distFaceY
    double* getVolStaggX() const;   // Returns volStaggX
    double* getVolStaggY() const;   // Returns volStaggY

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
    double atDistFaceX(int) const;      // Returns distFaceX[i] unsafely (does not check if 0 <= i < nx)
    double atDistFaceY(int) const;      // Returns distFaceY[j] unsafely (does not check if 0 <= j < ny)
    double atVolStaggX(int, int) const; // Returns volStaggX[i+j*(nx+1)] unsafely (does not check if 0 <= i < nx+1 and 0 <= j < ny+2)
    double atVolStaggY(int, int) const; // Returns volStaggY[i+j*(nx+2)] unsafely (does not check if 0 <= i < nx+2 and 0 <= j < ny+1)

    // Info functions
    void printBasicData() const;        // Prints: built, lx, ly, lz, nx, ny, nz=1
    void printNodeLocation() const;     // Prints node location and number (not nodeX nor nodeY)
    void printNodeDistances() const;    // Prints distX and distY
    void printSurfaces() const;         // Prints surfX and surfY
    void printVolumes() const;          // Prints vol
    void printFaceDistances() const;    // Prints distFaceX and distFaceY
    void printStaggeredVolumes() const; // Prints volStaggX and volStaggY
    void printMeshData() const;         // Prints all the previous member variables

    void saveMeshData() const;  // Saves the mesh parameters to different files to plot them later on with gnuplot
};

#endif /* NCMESH_H */
