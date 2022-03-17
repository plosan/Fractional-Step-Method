#ifndef NCMESH_H
#define NCMESH_H

#include <iostream>

class NCMesh {

private:
    bool built;     // Whether the mesh is built (functions) or not. Used to tell if getters can access member variables
    double lx;      // Length of the domain in the X coordinate
    double ly;      // Length of the domain in the Y coordinate
    double lz;      // Length of the domain in the Z coordinate
    int nx;         // Number of control volumes for discretisation in the X coordinate
    int ny;         // Number of control volumes for discretisation in the Y coordinate
    double* nodeX;  // Position of nodes in the X coordinate. Size: nx+2
    double* nodeY;  // Position of nodes in the Y coordinate. Size: ny+2
    double* faceX;  // Position in the X coordinate of the faces perpendicular to the X axis. Size: nx+1
    double* faceY;  // Position in the Y coordinate of the faces perpendicular to the Y axis. Size: ny+1
    double* surfX;  // Surface of the faces perpendicular to the X axis. Size: ny
    double* surfY;  // Surface of the faces perpendicular to the Y axis. Size: nx
    double* vol;    // Volume of the control volumes. Size: nx*ny

    bool computeFaceXY();
    bool computeNodeXY();
    bool computeSurfXY();
    bool computeVol();

public:

    // Constructors
    NCMesh();
    NCMesh(double _lx, double _ly, double _lz, int _nx, int _ny);
    NCMesh(double _lx, double _ly, double _lz, int _nx, int _ny, double* _faceX, double* _faceY);

    // Info functions
    void printMeshData() const; // Prints the mesh parameters
    void saveMeshData() const;  // Saves the mesh parameters to different files to plot it later on



};

#endif /* NCMESH_H */
