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
    double* distX;  // Distance between nodes in X coordinate. Size: nx-1
    double* distY;  // Distance between nodes in Y coordinate. Size: ny-1
    double* surfX;  // Surface of the faces perpendicular to the X axis. Size: ny
    double* surfY;  // Surface of the faces perpendicular to the Y axis. Size: nx
    double* vol;    // Volume of each control volume. Size: nx*ny

    // Compute grid parameters
    void computeNodeXY();   // Compute the position of nodes in X and Y coordinates
    void computeDistXY();   // Compute the distances between nodes in X coordinate and Y coordinate
    void computeFaceXY();   // Compute the location of faces perpendicular to the X axis and Y axis
    void computeSurfXY();   // Compute the surface of the faces perpendicular to the X axis and Y axis
    void computeVol();      // Compute the volume of each control volume


public:
    // Constructors
    RCGrid();
    RCGrid(double _lx, double _ly, double _lz, int _nx, int _ny);
    RCGrid(double _lx, double _ly, double _lz, int _nx, int _ny, double* _nodeX, double* _nodeY);

    // Getters
    int getNX() const;
    int getNY() const;
    double* getNodeX() const;
    double* getNodeY() const;
    double* getFaceX() const;
    double* getFaceY() const;
    


    // Info functions
    void printMeshData() const; // Prints the mesh parameters
    void saveMeshData() const;  // Saves the mesh parameters to different files to plot it later on

};
