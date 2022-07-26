#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cfloat>
#include <cstring>
#include <chrono>

#include "fsm.h"
#include "schemes.h"
#include "solver.h"


void allocateFluidVariables(const int nx, const int ny, double* &u, double* &v, double* &p) {
    /*
    allocateFluidVariables: allocates the arrays for x-velocity, y-velocity and pressure
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx    Control volumes in x direction          [const int]
        - ny    Control volumes in y direction          [const int]
        - u     X-velocity array. Size: (nx+1)*(ny+2)   [double* &]
        - v     Y-velocity array. Size: (nx+2)*(ny+1)   [double* &]
        - p     Pressure. Size: (nx+2)*(ny+2)           [double* v]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - none
    */
    // X-component of velocity
    u = (double*) calloc((nx+1)*(ny+2), sizeof(double));
    if(!u) {
        printf("Error: could not allocate enough memory for u\n");
        return;
    }
    // Y-component of velocity
    v = (double*) calloc((nx+2)*(ny+1), sizeof(double));
    if(!v) {
        printf("Error: could not allocate enough memory for v\n");
        return;
    }
    // Pressure
    p = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    if(!p) {
        printf("Error: could not allocate enough memory for p\n");
        return;
    }
}

void allocateOperatorR(const int nx, const int ny, double* &Ru, double* &Rv, double* &Ru_prev, double* &Rv_prev) {
    /*
    allocateOperatorR: allocates the arrays for the operator R, namely Ru, Rv, Ru_prev and Rv_prev
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx        Control volumes in x direction                                          [const int]
        - ny        Control volumes in y direction                                          [const int]
        - Ru        X-component of operator R at the current time. Size: (nx+1)*(ny+2)      [double* &]
        - Rv        Y-component of operator R at the current time. Size: (nx+2)*(ny+1)      [double* &]
        - Ru_prev   X-component of operator R at the previous time. Size: (nx+1)*(ny+2)     [double* &]
        - Rv_prev   Y-component of operator R at the previous time. Size: (nx+2)*(ny+1)     [double* &]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - none
    */
    // Operator Ru at time n
    Ru = (double*) calloc((nx+1)*(ny+2), sizeof(double));
    if(!Ru) {
        printf("Error: could not allocate enough memory for Ru\n");
        return;
    }

    // Operator Rv at time n
    Rv = (double*) calloc((nx+2)*(ny+1), sizeof(double));
    if(!Rv) {
        printf("Error: could not allocate enough memory for Rv\n");
        return;
    }

    // Operator Ru at time n-1
    Ru_prev = (double*) calloc((nx+1)*(ny+2), sizeof(double));
    if(!Ru_prev) {
        printf("Error: could not allocate enough memory for Ru\n");
        return;
    }

    // Operator Rv at time n-1
    Rv_prev = (double*) calloc((nx+2)*(ny+1), sizeof(double));
    if(!Rv) {
        printf("Error: could not allocate enough memory for Rv_prev\n");
        return;
    }
}


void allocatePredictorVelocities(const int nx, const int ny, double* &u_pred, double* &v_pred) {
    // X-velocity predictor
    u_pred = (double*) calloc((nx+1)*(ny+2), sizeof(double));
    if(!u_pred) {
        printf("Error: could not allocate enough memory for u_pred\n");
        return;
    }
    // Y-velocity predictor
    v_pred = (double*) calloc((nx+2)*(ny+1), sizeof(double));
    if(!v_pred) {
        printf("Error: could not allocate enough memory for u_pred\n");
        return;
    }
}


void allocateLinearSystemVariables(const int nx, const int ny, double* &A, double* &b) {
    // Linear system matrix
    A = (double*) calloc(5*(nx+2)*(ny+2), sizeof(double));
    if(!A) {
        printf("Error: could not allocate enough memory for coefficients matrix\n");
        return;
    }
    // Linear system vector
    b = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    if(!b) {
        printf("Error: could not allocate enough memory for independent terms vector\n");
        return;
    }
}


// Mass flows at faces
void computeMassFlowsStaggX(double* mx, double* my, const NCMesh m, const double* u, const double* v, const Properties props) {

    // Mesh sizes
    int nx = m.getNX(); // X-axis control volume count
    int ny = m.getNY(); // Y-axis control volume count

    // Mass flow through faces perpendicular to the X axis
    // For j = 0 and j = ny+1, mass flows are zero
    // mx = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    for(int j = 1; j < ny+1; j++) {
        double Ax = m.atSurfX(j);
        mx[j*(nx+2)] = props.rho * u[j*(nx+1)] * Ax;
        mx[j*(nx+2)+(nx+1)] = props.rho * u[j*(nx+1)+nx] * Ax;
        for(int i = 1; i < nx+1; i++) {
            double uW = u[j*(nx+1)+i-1];
            double uP = u[j*(nx+1)+i];
            mx[j*(nx+2)+i] = 0.5 * props.rho * (uP + uW) * Ax;
        }
    }

    // Mass flow through the faces perpendicular to the Y axis
    // For i = 0 and i = nx, mass flows are not needed
    // my = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    for(int i = 1; i < nx; i++) {
        double Ay_left = m.atSemiSurfY(i,1);        // Y-area to the left of the wall
        double Ay_right = m.atSemiSurfY(i+1,0);     // Y-area to the right of the wall
        for(int j = 0; j < ny+1; j++) {
            double vy_left = v[j*(nx+2)+i];         // Y-velocity to the left of the wall
            double vy_right = v[j*(nx+2)+i+1];      // Y-velocity to the right of the wall
            my[j*(nx+1)+i] = props.rho * (vy_left * Ay_left + vy_right * Ay_right);
        }
    }

}

void computeMassFlowsStaggY(double* mx, double* my, const NCMesh m, const double* u, const double* v, const Properties props) {

    // Mesh sizes
    int nx = m.getNX(); // X-axis control volume count
    int ny = m.getNY(); // Y-axis control volume count

    // Mass flow through faces perpendicular to the X axis
    // For j = 0 and j = ny, mass flows are not needed
    // mx = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    for(int j = 1; j < ny; j++) {
        double Ax_up = m.atSemiSurfX(j+1,0);    // X-area above the wall
        double Ax_down = m.atSemiSurfX(j,1);    // X-area below the wall
        for(int i = 0; i < nx+1; i++) {
            double ux_up = u[(j+1)*(nx+1)+i];   // X-velocity above the wall
            double ux_down = u[j*(nx+1)+i];     // X-velocity below the wall
            mx[j*(nx+1)+i] = props.rho * (ux_up * Ax_up + ux_down * Ax_down);
        }
    }

    // Mass flow through faces perpendicular to the Y axis
    // For j = 0 and j = ny+1, mass flows are zero
    // my = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    for(int i = 1; i < nx+1; i++) {
        double Ay = m.atSurfY(i);
        my[i] = props.rho * v[i] * Ay;
        my[(ny+1)*(nx+2)+i] = props.rho * v[(ny+1)*(nx+2)+i] * Ay;
        for(int j = 1; j < ny+1; j++) {
            double vS = v[(j-1)*(nx+2)+i];
            double vP = v[j*(nx+2)+i];
            my[j*(nx+2)+i] = 0.5 * props.rho * (vP + vS) * Ay;
        }
    }
}

void computeVelocitiesStaggX_CDS(double* ue, double* un, const int nx, const int ny, const double* u) {

    // Velocities uw and ue at the X-staggered control volume faces (west and east faces)
    // For j = 0 and j = ny+1 these velocities are not required
    // For i = 0 and i = nx+1 these velocities are equal to the velocity at the node
    // ue = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    for(int j = 1; j < ny+1; j++) {
        ue[j*(nx+2)] = u[j*(nx+1)];             // When i = 0
        ue[j*(nx+2)+(nx+1)] = u[j*(nx+1)+nx];   // When i = nx + 1
        for(int i = 1; i < nx+1; i++) {
            double u_left = u[j*(nx+1)+i-1];
            double u_right = u[j*(nx+1)+i];
            ue[j*(nx+2)+i] = 0.5 * (u_left + u_right);
        }
    }

    // Velocities us and un at the X-staggered control volumes faces (south and north faces)
    // For i = 0 and i = nx these velocities are not required
    // For j = 0 or j = ny and 1 <= i < nx, these velocities are equal to the velocity at the node
    // un = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    for(int i = 1; i < nx; i++) {
        un[i] = u[i];                       // When j = 0
        un[ny*(nx+1)+i] = u[(ny+1)*(nx+1)+i];   // When j = ny
        for(int j = 1; j < ny; j++) {
            double u_below = u[j*(nx+1)+i];
            double u_above = u[(j+1)*(nx+1)+i];
            un[j*(nx+1)+i] = 0.5 * (u_below + u_above);
        }
    }
}

void computeVelocitiesStaggY_CDS(double* ve, double* vn, const int nx, const int ny, const double* v) {

    // Velocities vw and ve at the Y-staggered control volume faces (west and east)
    // For j = 0 and j = ny these velocities are not required
    // For 1 <= j < ny and (i = 0 or i = nx), these velocities are equal to those at the nodes
    // ve = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    for(int j = 1; j < ny; j++) {
        ve[j*(nx+1)] = v[j*(nx+2)];             // When i = 0
        ve[j*(nx+1)+nx] = v[j*(nx+2)+(nx+1)];   // When i = nx
        for(int i = 1; i < nx; i++) {
            double v_left = v[j*(nx+2)+i];
            double v_right = v[j*(nx+2)+i+1];
            ve[j*(nx+1)+i] = 0.5 * (v_left + v_right);
        }
    }

    // Velocities vs and vn at the Y-staggered control volume faces (south and north)
    // For i = 0 and i = nx+1 these velocities are not required
    // For j = 0 or j = ny+1 and 1 <= i < nx+1, these velocities are equal to those at the nodes
    // vn = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    for(int i = 1; i < nx+1; i++) {
        vn[i] = v[i];                               // When j = 0
        vn[(ny+1)*(nx+2)+i] = v[(ny+1)*(nx+2)+i];   // When j = ny+1
        for(int j = 1; j < ny+1; j++) {
            double v_below = v[(j-1)*(nx+2)+i];
            double v_above = v[j*(nx+2)+i];
            vn[j*(nx+2)+i] = 0.5 * (v_below + v_above);
        }
    }
}

void computeVelocitiesStaggX_QUICK(double* ue, double* un, const NCMesh m, const double* mx, const double* u) {

    // Mesh sizes
    int nx = m.getNX(); // X-axis control volume count
    int ny = m.getNY(); // Y-axis control volume count

    // Velocities uw and ue at the X-staggered control volume faces (west and east faces)
    // For j = 0 and j = ny+1 these velocities are not required
    // For i = 0 and i = nx+1 these velocities are equal to the velocity at the node
    // ue = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    for(int j = 1; j < ny+1; j++) {
        ue[j*(nx+2)] = u[j*(nx+1)];             // When i = 0
        ue[j*(nx+2)+(nx+1)] = u[j*(nx+1)+nx];   // When i = nx + 1
        for(int i = 1; i < nx+1; i++) {

            double mf = mx[j*(nx+2)+i];         // Mass flow
            if(mf > 0) {    // From left to right
                if(i == 1) {    // There is no second-upwind node (in this case W), so we use CDS instead
                    double u_left = u[j*(nx+1)+i-1];
                    double u_right = u[j*(nx+1)+i];
                    ue[j*(nx+2)+i] = 0.5 * (u_left + u_right);
                } else {
                    // Property values
                    double u_UU = u[j*(nx+1)+i-2];  // X-velocity at the second upstream node
                    double u_U = u[j*(nx+1)+i-1];   // X-velocity at the first upstream node
                    double u_D = u[j*(nx+1)+i];     // X-velocity at the first downstream node
                    // Nodes and face positions
                    double x_UU = m.atFaceX(i-2);   // X-position of the second upstream face
                    double x_U = m.atFaceX(i-1);    // X-position of the first upstream face
                    double x_D = m.atFaceX(i);      // X-position of the first downstream face
                    double xe = m.atNodeX(i);       // X-position of the node where the convective property is being computed
                    // QUICK scheme
                    ue[j*(nx+2)+i] = schemeQUICK(u_D, u_U, u_UU, x_D, x_U, x_UU, xe);
                }
            } else if(mf < 0) {
                if(i == nx) {   // There is not second-upwind node (in this case EE), so we use CDS instead
                    double u_left = u[j*(nx+1)+i-1];
                    double u_right = u[j*(nx+1)+i];
                    ue[j*(nx+2)+i] = 0.5 * (u_left + u_right);
                } else {
                    // Property values
                    double u_UU = u[j*(nx+1)+i+1];  // X-velocity at the second upstream node
                    double u_U = u[j*(nx+1)+i];     // X-velocity at the first upstream node
                    double u_D = u[j*(nx+1)+i-1];   // X-velocity at the first downstream node
                    // Nodes and face positions
                    double x_UU = m.atFaceX(i+1);   // X-position of the second upstream face
                    double x_U = m.atFaceX(i);      // X-position of the first upstream face
                    double x_D = m.atFaceX(i-1);    // X-position of the first downstream face
                    double xe = m.atNodeX(i);       // X-position of the node where the convective property is being computed
                    // QUICK scheme
                    ue[j*(nx+2)+i] = schemeQUICK(u_D, u_U, u_UU, x_D, x_U, x_UU, xe);
                }
            } else {    // Rare case, then use CDS
                double u_left = u[j*(nx+1)+i-1];
                double u_right = u[j*(nx+1)+i];
                ue[j*(nx+2)+i] = 0.5 * (u_left + u_right);
            }

        }
    }

    // Velocities us and un at the X-staggered control volumes faces (south and north faces)
    // For i = 0 and i = nx these velocities are not required
    // For j = 0 or j = ny and 1 <= i < nx, these velocities are equal to the velocity at the node
    // un = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    for(int i = 1; i < nx; i++) {
        un[i] = u[i];                       // When j = 0
        un[ny*(nx+1)+i] = u[(ny+1)*(nx+1)+i];   // When j = ny
        for(int j = 1; j < ny; j++) {
            double u_below = u[j*(nx+1)+i];
            double u_above = u[(j+1)*(nx+1)+i];
            un[j*(nx+1)+i] = 0.5 * (u_below + u_above);
        }
    }
}

void computeVelocitiesStaggY_QUICK(double* ve, double* vn, const NCMesh m, const double* my, const double* v) {

    // Mesh sizes
    int nx = m.getNX(); // X-axis control volume count
    int ny = m.getNY(); // Y-axis control volume count

    // Velocities vw and ve at the Y-staggered control volume faces (west and east)
    // For j = 0 and j = ny these velocities are not required
    // For 1 <= j < ny and (i = 0 or i = nx), these velocities are equal to those at the nodes
    // ve = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    for(int j = 1; j < ny; j++) {
        ve[j*(nx+1)] = v[j*(nx+2)];             // When i = 0
        ve[j*(nx+1)+nx] = v[j*(nx+2)+(nx+1)];   // When i = nx
        for(int i = 1; i < nx; i++) {
            double v_left = v[j*(nx+2)+i];
            double v_right = v[j*(nx+2)+i+1];
            ve[j*(nx+1)+i] = 0.5 * (v_left + v_right);
        }
    }

    // Velocities vs and vn at the Y-staggered control volume faces (south and north)
    // For i = 0 and i = nx+1 these velocities are not required
    // For j = 0 or j = ny+1 and 1 <= i < nx+1, these velocities are equal to those at the nodes
    // vn = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    for(int i = 1; i < nx+1; i++) {
        vn[i] = v[i];                               // When j = 0
        vn[(ny+1)*(nx+2)+i] = v[(ny+1)*(nx+2)+i];   // When j = ny+1
        for(int j = 1; j < ny+1; j++) {

            double mf = my[j*(nx+2)+i];      // Mass flow at face

            if(mf > 0) {
                if(j == 1) {    // There is no second upwind node (in this case S), so we use CDS instead
                    double v_below = v[(j-1)*(nx+2)+i];
                    double v_above = v[j*(nx+2)+i];
                    vn[j*(nx+2)+i] = 0.5 * (v_below + v_above);
                } else {
                    // Property values
                    double v_UU = v[(j-2)*(nx+2)+i];    // Y-velocity at the second upstream node
                    double v_U = v[(j-1)*(nx+2)+i];     // Y-velocity at the first upstream node
                    double v_D = v[j*(nx+2)+i];         // Y-velocity at the first downstream node
                    // Nodes and face positions
                    double y_UU = m.atFaceY(j-2);       // Y-position of the second upstream face
                    double y_U = m.atFaceY(j-1);        // Y-position of the first upstream face
                    double y_D = m.atFaceY(j);          // Y-position of the first downstream face
                    double yn = m.atNodeY(j);           // Y-position of the node where the convective property is being computed
                    // QUICK scheme
                    vn[j*(nx+2)+i] = schemeQUICK(v_D, v_U, v_UU, y_D, y_U, y_UU, yn);
                }
            } else if(mf < 0) {
                if(j == ny) {   // There is no second upwind node (in this case NN), so we use CDS instead
                    double v_below = v[(j-1)*(nx+2)+i];
                    double v_above = v[j*(nx+2)+i];
                    vn[j*(nx+2)+i] = 0.5 * (v_below + v_above);
                } else {
                    // Property values
                    double v_UU = v[(j+1)*(nx+2)+i];    // Y-velocity at the second upstream node
                    double v_U = v[j*(nx+2)+i];         // Y-velocity at the first upstream node
                    double v_D = v[(j-1)*(nx+2)+i];     // Y-velocity at the first downstream node
                    // Nodes and face positions
                    double y_UU = m.atFaceY(j+1);       // Y-position of the second upstream face
                    double y_U = m.atFaceY(j);          // Y-position of the first upstream face
                    double y_D = m.atFaceY(j-1);        // Y-position of the first downstream face
                    double yn = m.atNodeY(j);           // Y-position of the node where the convective property is being computed
                    // QUICK scheme
                    vn[j*(nx+2)+i] = schemeQUICK(v_D, v_U, v_UU, y_D, y_U, y_UU, yn);
                }
            } else {    // Rare case, then use CDS
                double v_below = v[(j-1)*(nx+2)+i];
                double v_above = v[j*(nx+2)+i];
                vn[j*(nx+2)+i] = 0.5 * (v_below + v_above);
            }

        }
    }
}

void computeRu2(double* Ru, const NCMesh m, const double* u, const double* v, const Properties props) {

    // Mesh sizes
    int nx = m.getNX(); // X-axis control volume count
    int ny = m.getNY(); // Y-axis control volume count

    // Mass flows at faces for the X-staggered mesh
    double* mx = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    double* my = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    computeMassFlowsStaggX(mx, my, m, u, v, props);

    // Velocities at faces for the X-staggered mesh
    double* u_hor = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    double* u_ver = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    computeVelocitiesStaggX_CDS(u_hor, u_ver, nx, ny, u);

    for(int j = 1; j < ny+1; j++) {
        for(int i = 1; i < nx; i++) {
            // Nodal velocities
            int node = j * (nx + 1) + i;
            double uP = u[node];
            double uW = u[node-1];
            double uE = u[node+1];
            double uS = u[node-(nx+1)];
            double uN = u[node+(nx+1)];
            // Faces velocities
            double uw = u_hor[j*(nx+2)+i];
            double ue = u_hor[j*(nx+2)+i+1];
            double us = u_ver[(j-1)*(nx+1)+i];
            double un = u_ver[j*(nx+1)+i];
            // Areas
            double Ax = m.atSurfX(j);
            double Ay = m.atSurfY_StaggX(i);
            // Mass flows
            double mw = mx[j*(nx+2)+i];
            double me = mx[j*(nx+2)+i+1];
            double mn = my[j*(nx+1)+i];
            double ms = my[(j-1)*(nx+1)+i];
            // Operator R(u)
            double integral1 = -(me * ue - mw * uw + mn * un - ms * us);
            double integral2 = Ax * (uE - uP) / m.atDistFaceX(i) - Ax * (uP - uW) / m.atDistFaceX(i-1);
            integral2 += Ay * (uN - uP) / m.atDistY(j) - Ay * (uP - uS) / m.atDistY(j-1);
            integral2 *= props.mu;
            Ru[j*(nx+1)+i] = (integral1 + integral2) / m.atVolStaggX(i,j);
        }
    }

    // Free allocated memory
    free(mx);
    free(my);
    free(u_hor);
    free(u_ver);
}

void computeRu3(double* Ru, const NCMesh m, const double* u, const double* ux_staggX, const double* uy_staggX, const double* mx_staggX, const double* my_staggX, const Properties props) {

    // Mesh sizes
    int nx = m.getNX(); // X-axis control volume count
    int ny = m.getNY(); // Y-axis control volume count

    for(int j = 1; j < ny+1; j++) {
        for(int i = 1; i < nx; i++) {
            // Nodal velocities
            int node = j * (nx + 1) + i;
            double uP = u[node];
            double uW = u[node-1];
            double uE = u[node+1];
            double uS = u[node-(nx+1)];
            double uN = u[node+(nx+1)];
            // Faces velocities
            double uw = ux_staggX[j*(nx+2)+i];
            double ue = ux_staggX[j*(nx+2)+i+1];
            double us = uy_staggX[(j-1)*(nx+1)+i];
            double un = uy_staggX[j*(nx+1)+i];
            // Areas
            double Ax = m.atSurfX(j);
            double Ay = m.atSurfY_StaggX(i);
            // Mass flows
            double mw = mx_staggX[j*(nx+2)+i];
            double me = mx_staggX[j*(nx+2)+i+1];
            double mn = my_staggX[j*(nx+1)+i];
            double ms = my_staggX[(j-1)*(nx+1)+i];
            // Operator R(u)
            double integral1 = -(me * ue - mw * uw + mn * un - ms * us);
            double integral2 = Ax * (uE - uP) / m.atDistFaceX(i) - Ax * (uP - uW) / m.atDistFaceX(i-1);
            integral2 += Ay * (uN - uP) / m.atDistY(j) - Ay * (uP - uS) / m.atDistY(j-1);
            integral2 *= props.mu;
            Ru[j*(nx+1)+i] = (integral1 + integral2) / m.atVolStaggX(i,j);
        }
    }

}

void computeRv3(double* Rv, const NCMesh m, const double* v, const double* vx_staggY, const double* vy_staggY, const double* mx_staggY, const double* my_staggY, const Properties props) {

    // Mesh sizes
    int nx = m.getNX(); // X-axis control volume count
    int ny = m.getNY(); // Y-axis control volume count

    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny; j++) {
            // Nodal velocities
            int node = j * (nx + 2) + i;
            double vP = v[node];
            double vW = v[node-1];
            double vE = v[node+1];
            double vS = v[node-(nx+2)];
            double vN = v[node+(nx+2)];
            // Face velocities
            double vw = vx_staggY[j*(nx+1)+i-1];
            double ve = vx_staggY[j*(nx+1)+i];
            double vs = vy_staggY[j*(nx+2)+i];
            double vn = vy_staggY[(j+1)*(nx+2)+i];
            // Areas
            double Ax = m.atSurfX_StaggY(j);
            double Ay = m.atSurfY(i);
            // Mass flows
            double mw = mx_staggY[j*(nx+1)+i-1];
            double me = mx_staggY[j*(nx+1)+i];
            double ms = my_staggY[j*(nx+2)+i];
            double mn = my_staggY[(j+1)*(nx+2)+i];
            // Operator R(v)
            double integral1 = -(me * ve - mw * vw + mn * vn - ms * vs);
            double integral2 = Ax * (vE - vP) / m.atDistX(i) - Ax * (vP - vW) / m.atDistX(i-1);
            integral2 += Ay * (vN - vP) / m.atDistFaceY(j) - Ay * (vP - vS) / m.atDistFaceY(j-1);
            integral2 *= props.mu;
            Rv[j*(nx+2)+i] = (integral1 + integral2) / m.atVolStaggY(i,j);
        }
    }



}

void computeRv2(double* Rv, const NCMesh m, const double* u, const double* v, const Properties props) {

    // Mesh sizes
    int nx = m.getNX(); // X-axis control volume count
    int ny = m.getNY(); // Y-axis control volume count

    // Mass flows at faces for the Y-staggered mesh
    double* mx = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    double* my = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    computeMassFlowsStaggY(mx, my, m, u, v, props);

    // Velocities at faces for the Y-staggered mesh
    double* v_hor = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    double* v_ver = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    computeVelocitiesStaggY_CDS(v_hor, v_ver, nx, ny, v);

    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny; j++) {
            // Nodal velocities
            int node = j * (nx + 2) + i;
            double vP = v[node];
            double vW = v[node-1];
            double vE = v[node+1];
            double vS = v[node-(nx+2)];
            double vN = v[node+(nx+2)];
            // Face velocities
            double vw = v_hor[j*(nx+1)+i-1];
            double ve = v_hor[j*(nx+1)+i];
            double vs = v_ver[j*(nx+2)+i];
            double vn = v_ver[(j+1)*(nx+2)+i];
            // Areas
            double Ax = m.atSurfX_StaggY(j);
            double Ay = m.atSurfY(i);
            // Mass flows
            double mw = mx[j*(nx+1)+i-1];
            double me = mx[j*(nx+1)+i];
            double ms = my[j*(nx+2)+i];
            double mn = my[(j+1)*(nx+2)+i];
            // Operator R(v)
            double integral1 = -(me * ve - mw * vw + mn * vn - ms * vs);
            double integral2 = Ax * (vE - vP) / m.atDistX(i) - Ax * (vP - vW) / m.atDistX(i-1);
            integral2 += Ay * (vN - vP) / m.atDistFaceY(j) - Ay * (vP - vS) / m.atDistFaceY(j-1);
            integral2 *= props.mu;
            Rv[j*(nx+2)+i] = (integral1 + integral2) / m.atVolStaggY(i,j);
        }
    }

    free(mx);
    free(my);
    free(v_hor);
    free(v_ver);
}

void computeRu(double* Ru, const NCMesh m, const double* u, const double* v, const Properties props) {

    int nx = m.getNX();
    int ny = m.getNY();

    for(int j = 1; j < ny+1; j++) {
        for(int i = 1; i < nx; i++) {
            // Nodal velocities
            int node = j * (nx + 1) + i;
            double uP = u[node];
            double uW = u[node-1];
            double uE = u[node+1];
            double uS = u[node-(nx+1)];
            double uN = u[node+(nx+1)];
            // Faces velocities
            double uw = schemeCDS(uP, uW);
            double ue = schemeCDS(uP, uE);
            double us = (j > 1 ? 0.5*(uP + uS) : uS);
            double un = (j < ny ? 0.5*(uP + uN) : uN);
            // Areas
            double Ax = m.atSurfX(j);
            double Ay = m.atSurfY_StaggX(i);
            double Ay_left = m.atSemiSurfY(i,1);
            double Ay_right = m.atSemiSurfY(i+1,0);
            // Mass flows
            double mw = 0.5 * props.rho * (uP + uW) * Ax;
            double me = 0.5 * props.rho * (uP + uE) * Ax;
            double mn = props.rho * (v[j*(nx+2)+i] * Ay_left + v[j*(nx+2)+i+1] * Ay_right);
            double ms = props.rho * (v[(j-1)*(nx+2)+i] * Ay_left + v[(j-1)*(nx+2)+i+1] * Ay_right);
            // Operator R(u)
            double integral1 = -(me * ue - mw * uw + mn * un - ms * us);
            double integral2 = Ax * (uE - uP) / m.atDistFaceX(i) - Ax * (uP - uW) / m.atDistFaceX(i-1);
            integral2 += Ay * (uN - uP) / m.atDistY(j) - Ay * (uP - uS) / m.atDistY(j-1);
            integral2 *= props.mu;
            Ru[j*(nx+1)+i] = (integral1 + integral2) / m.atVolStaggX(i,j);
        }
    }
}

void computeRv(double* Rv, const NCMesh m, const double* u, const double* v, const Properties props) {

    int nx = m.getNX();
    int ny = m.getNY();

    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny; j++) {
            // Nodal velocities
            int node = j * (nx + 2) + i;
            double vP = v[node];
            double vW = v[node-1];
            double vE = v[node+1];
            double vS = v[node-(nx+2)];
            double vN = v[node+(nx+2)];
            // Face velocities
            double vw = (i > 1 ? 0.5*(vP + vW) : vW);
            double ve = (i < nx ? 0.5*(vP + vE) : vE);
            double vs = schemeCDS(vP, vS);
            double vn = schemeCDS(vP, vN);
            // Areas
            double Ax = m.atSurfX_StaggY(j);
            double Ax_up = m.atSemiSurfX(j+1,0);
            double Ax_down = m.atSemiSurfX(j,1);
            double Ay = m.atSurfY(i);
            // Mass flows
            double mw = props.rho * (u[j*(nx+1)+i-1] * Ax_down + u[(j+1)*(nx+1)+i-1] * Ax_up);
            double me = props.rho * (u[j*(nx+1)+i] * Ax_down + u[(j+1)*(nx+1)+i] * Ax_up);
            double ms = 0.5 * props.rho * (vP + vS) * Ay;
            double mn = 0.5 * props.rho * (vP + vN) * Ay;
            // Operator R(v)
            double integral1 = -(me * ve - mw * vw + mn * vn - ms * vs);
            double integral2 = Ax * (vE - vP) / m.atDistX(i) - Ax * (vP - vW) / m.atDistX(i-1);
            integral2 += Ay * (vN - vP) / m.atDistFaceY(j) - Ay * (vP - vS) / m.atDistFaceY(j-1);
            integral2 *= props.mu;
            Rv[j*(nx+2)+i] = (integral1 + integral2) / m.atVolStaggY(i,j);
        }
    }
}

void computePredictorVelocityU(double* u_pred, const int nx, const int ny, const double* u, const double* Ru, const double* Ru_prev, const Properties props, const double tstep) {
    for(int i = 1; i < nx; i++) {
        for(int j = 1; j < ny+1; j++) {
            int node = j * (nx + 1) + i;
            u_pred[node] = u[node] + (1.5 * Ru[node] - 0.5 * Ru_prev[node]) * tstep / props.rho;
        }
    }
}


void computePredictorVelocityV(double* v_pred, const int nx, const int ny, const double* v, const double* Rv, const double* Rv_prev, const Properties props, const double tstep) {
    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny; j++) {
            int node = j * (nx + 2) + i;
            v_pred[node] = v[node] + (1.5 * Rv[node] - 0.5 * Rv_prev[node]) * tstep / props.rho;
        }
    }
}

void computeDiscretizationCoefficients(double* A, double* b, const NCMesh m, const double* u_pred, const double* v_pred, const Properties props, const double tstep) {

    // Mesh sizes
    int nx = m.getNX(); // X-axis control volume count
    int ny = m.getNY(); // Y-axis control volume count

    // Compute discretization coefficients
    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny+1; j++) {
            // Surfaces and distances
            double Ax = m.atSurfX(j);       // X-surface
            double Ay = m.atSurfY(i);       // Y-surface
            double dPE = m.atDistX(i);      // Current node - East node distance
            double dPW = m.atDistX(i-1);    // Current node - West node distance
            double dPN = m.atDistY(j);      // Current node - North node distance
            double dPS = m.atDistY(j-1);    // Current node - South node distance
            // Discretization coefficients
            int node = j * (nx + 2) + i;    // Node number
            A[5*node] = Ay / dPS;           // South node
            A[5*node+1] = Ax / dPW;         // West node
            A[5*node+2] = Ax / dPE;         // East node
            A[5*node+3] = Ay / dPN;         // North node
            A[5*node+4] = A[5*node] + A[5*node+1] + A[5*node+2] + A[5*node+3];  // Current node
            // Independent term
            int node_east = j * (nx + 1) + i;           // X-staggered east node number
            int node_west = j * (nx + 1) + i - 1;       // X-staggered west node number
            int node_north = j * (nx + 2) + i;          // Y-staggered north node number
            int node_south = (j - 1) * (nx + 2) + i;    // Y-staggered south node number
            b[node] = (-props.rho / tstep) * (u_pred[node_east] * Ax - u_pred[node_west] * Ax + v_pred[node_north] * Ay - v_pred[node_south] * Ay);    // Independent term
        }
    }
}

void computeNextVelocityField(double* u, double* v, const NCMesh m, const double* u_pred, const double* v_pred, const double* p, const Properties props, const double tstep, double& maxDiff) {

    // Mesh size
    int nx = m.getNX(); // Number of control volumes along x axis
    int ny = m.getNY(); // Number of control volumes along y axis

    // X-component of velocity
    for(int i = 1; i < nx; i++) {
        for(int j = 1; j < ny+1; j++) {
            int node = j * (nx + 2) + i;                                // Pressure node number
            double px = (p[node+1] - p[node]) / m.atDistX(i);           // Partial derivative of pressure with respect to x
            node = j * (nx + 1) + i;                                    // X-velocity node number
            double u_next = u_pred[node] - (tstep / props.rho) * px;    // X-velocity at time n+1
            maxDiff = std::max(maxDiff, std::abs(u[node] - u_next));    // Update maximum difference
            u[node] = u_next;                                           // Update velocity
        }
    }

    // Y-component of velocity
    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny; j++) {
            int node = j * (nx + 2) + i;                                // Pressure node number
            double py = (p[node+(nx+2)] - p[node]) / m.atDistY(j);      // Partial derivative of pressure with respect to y
            node = j * (nx + 2) + i;                                    // Y-velocity node number
            double v_next = v_pred[node] - (tstep / props.rho) * py;    // Y-velocity at time instant n+1
            maxDiff = std::max(maxDiff, std::abs(v[node] - v_next));    // Update maximum difference to check convergence
            v[node] = v_next;                                           // Update velocity
        }
    }
}

// void computeVelocityU(double* u, const NCMesh m, const double* u_pred, const double* p, const Properties props, const double tstep, double& maxDiff) {
//
//     // Mesh size
//     int nx = m.getNX(); // Number of control volumes along x axis
//     int ny = m.getNY(); // Number of control volumes along y axis
//
//     // X-component of velocity
//     for(int i = 1; i < nx; i++) {
//         for(int j = 1; j < ny+1; j++) {
//             int node = j * (nx + 2) + i;                                // Pressure node number
//             double px = (p[node+1] - p[node]) / m.atDistX(i);           // Partial derivative of pressure with respect to x
//             node = j * (nx + 1) + i;                                    // X-velocity node number
//             double u_next = u_pred[node] - (tstep / props.rho) * px;    // X-velocity at time n+1
//             maxDiff = std::max(maxDiff, std::abs(u[node] - u_next));    // Update maximum difference
//             u[node] = u_next;                                           // Update velocity
//         }
//     }
//
// }
//
// void computeVelocityV(double* v, const NCMesh m, const double* v_pred, const double* p, const Properties props, const double tstep, double& maxDiff) {
//
//     // Mesh size
//     int nx = m.getNX(); // Number of control volumes along x axis
//     int ny = m.getNY(); // Number of control volumes along y axis
//
//     // Y-component of velocity
//     for(int i = 1; i < nx+1; i++) {
//         for(int j = 1; j < ny; j++) {
//             int node = j * (nx + 2) + i;                                // Pressure node number
//             double py = (p[node+(nx+2)] - p[node]) / m.atDistY(j);      // Partial derivative of pressure with respect to y
//             node = j * (nx + 2) + i;                                    // Y-velocity node number
//             double v_next = v_pred[node] - (tstep / props.rho) * py;    // Y-velocity at time instant n+1
//             maxDiff = std::max(maxDiff, std::abs(v[node] - v_next));    // Update maximum difference to check convergence
//             v[node] = v_next;                                           // Update velocity
//         }
//     }
//
// }

void computeTimeStep(double &tstep, const NCMesh m, const double* u, const double* v, const Properties props, const double tol) {

    // Mesh size
    int nx = m.getNX(); // Number of control volumes along x axis
    int ny = m.getNY(); // Number of control volumes along y axis

    // Compute convective time limit
    double tc = DBL_MAX;    // Large initial value
    // Convective time limit for x axis
    for(int i = 0; i < nx+1; i++) {
        for(int j = 0; j < ny+2; j++) {
            double vel = std::abs(u[j*(nx+1)+i]);
            if(vel > tol)
                tc = std::min(tc, m.atDistX(i) / vel);
        }
    }
    // Convective time limit for y axis
    for(int i = 0; i < nx+2; i++) {
        for(int j = 0; j < ny+1; j++) {
            double vel = std::abs(v[j*(nx+2)+i]);
            if(vel > tol)
                tc = std::min(tc, m.atDistY(j) / vel);
        }
    }

    // Diffusive time limit
    double td = DBL_MAX;    // Large initial value
    double nu = props.mu / props.rho;   // Kinematic viscosity
    // Diffuse time limit for x axis
    for(int i = 0; i < nx+1; i++)
        td = std::min(td, m.atDistX(i) * m.atDistX(i) / nu);
    // Diffuse time limit for y axis
    for(int j = 0; j < ny+1; j++)
        td = std::min(td, m.atDistY(j) * m.atDistY(j) / nu);
    td *= 0.5;

    // Compute new time step
    tstep = 0.2 * std::min(tc, td);
}


void updateOperatorR(double* Ru_prev, double* Rv_prev, const double* Ru, const double* Rv, const int nx, const int ny) {
    // Update Ru
    for(int i = 1; i < nx; i++) {
        for(int j = 1; j < ny+1; j++) {
            int node = j * (nx + 1) + i;
            Ru_prev[node] = Ru[node];
        }
    }
    // Update Rv
    for(int i = 1; i < nx+1; i++) {
        for(int j = 1; j < ny; j++) {
            int node = j * (nx + 2) + i;
            Rv_prev[node] = Rv[node];
        }
    }
}

void computeCenteredNodesVelocities(double* u_col, double* v_col, const int nx, const int ny, const double* u, const double* v) {

    // COMPUTE U_COL
    // x-velocity at the center of the non-staggered control volumes
    // Left and right columns
    for(int j = 0; j < ny+2; j++) {
        u_col[j*(nx+2)] = u[j*(nx+1)];
        u_col[j*(nx+2)+nx+1] = u[j*(nx+1)+nx];
    }
    // Central columns
    for(int i = 1; i < nx+1; i++) {
        for(int j = 0; j < ny+2; j++) {
            double uW = u[j*(nx+1)+i-1];
            double uE = u[j*(nx+1)+i];
            u_col[j*(nx+2)+i] = schemeCDS(uW, uE);
        }
    }

    // Compute V_COL
    // y-velocity at the center of the non-staggered control volumes
    // Lower and upper rows
    for(int i = 0; i < nx+2; i++) {
        v_col[i] = v[i];
        v_col[(ny+1)*(nx+2)+i] = v[ny*(nx+2)+i];
    }
    // Central rows
    for(int j = 1; j < ny+1; j++) {
        for(int i = 0; i < nx+2; i++) {
            double vN = v[j*(nx+2)+i];
            double vS = v[(j-1)*(nx+2)+i];
            v_col[j*(nx+2)+i] = schemeCDS(vS, vN);
        }
    }


}


void printVariablesToFile(const NCMesh m, const double* u, const double* v, const double* p, const int Re, const int precision, const int time) {

    // Mesh size
    int nx = m.getNX(); // Number of control volumes along x axis
    int ny = m.getNY(); // Number of control volumes along y axis

    // X-velocities at the pressure nodes
    double* u_col = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    if(!u_col) {
        printf("Error: could not allocate enough memory for u_col\n");
        return;
    }

    // Y-velocities at the pressure nodes
    double* v_col = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    if(!v_col) {
        printf("Error: could not allocate enough memory for v_col\n");
        return;
    }

    // Compute velocities at the pressure nodes
    computeCenteredNodesVelocities(u_col, v_col, nx, ny, u, v);

    // Save velocity norm
    std::string filename = "../plots/vel_" + std::to_string(nx) + "_" + std::to_string(ny) + "_" + std::to_string(Re) + "_" + std::to_string(time) + ".txt";
    printVelocityToFile(m, u_col, v_col, filename.c_str(), 5);

    // Save x-velocity
    filename = "../plots/u_" + std::to_string(nx) + "_" + std::to_string(ny) + "_" + std::to_string(Re) + "_" + std::to_string(time) + ".txt";
    printVelocityUToFile(m, u_col, filename.c_str(), 5);

    // Save y-velocity
    filename = "../plots/v_" + std::to_string(nx) + "_" + std::to_string(ny) + "_" + std::to_string(Re) + "_" + std::to_string(time) + ".txt";
    printVelocityVToFile(m, v_col, filename.c_str(), 5);

    free(u_col);
    free(v_col);

    // Save pressure
    filename = "../plots/p_" + std::to_string(nx) + "_" + std::to_string(ny) + "_" + std::to_string(Re) + "_" + std::to_string(time) + ".txt";
    printPressureToFile(m, p, filename.c_str(), 5);
}

void printVelocityToFile(const NCMesh m, const double* u_col, const double* v_col, const char* filename, const int precision) {

    std::ofstream file;
    file.open(filename);

    if(!file.is_open()) {
        printf("Error. Could not open file\n");
        file.close();
        return;
    }

    int nx = m.getNX();
    int ny = m.getNY();

    file << std::setprecision(precision) << std::fixed;
    for(int i = 0; i < nx+2; i++) {
        for(int j = 0; j < ny+2; j++) {
            int node = j * (nx + 2) + i;
            double u = u_col[node];
            double v = v_col[node];
            double norm = std::sqrt(u * u + v * v);
            file << m.atNodeX(i) << " " << m.atNodeY(j) << " " << norm << " " << u << " " << v << std::endl;
        }
        file << std::endl;
    }


    file.close();

}

void printVelocityUToFile(const NCMesh m, const double* u_col, const char* filename, const int precision) {

    std::ofstream file;
    file.open(filename);

    if(!file.is_open()) {
        printf("Error. Could not open file\n");
        file.close();
        return;
    }

    int nx = m.getNX();
    int ny = m.getNY();

    file << std::setprecision(precision) << std::fixed;
    for(int i = 0; i < nx+2; i++) {
        for(int j = 0; j < ny+2; j++) {
            int node = j * (nx + 2) + i;
            file << m.atNodeX(i) << " " << m.atNodeY(j) << " " << u_col[node] << std::endl;
        }
        file << std::endl;
    }


    file.close();

}

void printVelocityVToFile(const NCMesh m, const double* v_col, const char* filename, const int precision) {

    std::ofstream file;
    file.open(filename);

    if(!file.is_open()) {
        printf("Error. Could not open file\n");
        file.close();
        return;
    }

    int nx = m.getNX();
    int ny = m.getNY();

    file << std::setprecision(precision) << std::fixed;
    for(int i = 0; i < nx+2; i++) {
        for(int j = 0; j < ny+2; j++) {
            int node = j * (nx + 2) + i;
            file << m.atNodeX(i) << " " << m.atNodeY(j) << " " << v_col[node] << std::endl;
        }
        file << std::endl;
    }


    file.close();
}


void printPressureToFile(const NCMesh m, const double* p, const char* filename, const int precision) {

    std::ofstream file;
    file.open(filename);

    if(!file.is_open()) {
        printf("Error. Could not open file\n");
        file.close();
        return;
    }

    int nx = m.getNX();
    int ny = m.getNY();

    file << std::setprecision(precision) << std::fixed;
    for(int i = 0; i < nx+2; i++) {
        for(int j = 0; j < ny+2; j++) {
            int node = j * (nx + 2) + i;
            file << m.atNodeX(i) << " " << m.atNodeY(j) << " " << p[node] << std::endl;
        }
        file << std::endl;
    }


    file.close();
}

// Lid-driven cavity

void lid_driven::mainLoop(const double rho, const double mu, const double u_ref, const double p_ref, const double L, const int nx, const int ny, const double tstep0, const int maxIt, const double tol, const double sstol) {

    // General variables
    Properties props = {rho, mu};   // Fluid properties
    NCMesh m(L, L, 1, nx, ny);      // Node centered staggered mesh
    double t = 0;           // Current time
    int it = 0;             // Iteration counter
    bool steady = false;    // Steady state bool. true = steady state reached, false = steady state not reached
    double tstep = tstep0;  // Initial time step

    // Fluid variables
    double* u;  // X-velocity array
    double* v;  // Y-velocity array
    double* p;  // Pressure array
    allocateFluidVariables(nx, ny, u, v, p);                // Allocate the previous arrays
    lid_driven::setInitialMaps(u, v, p, m, u_ref, p_ref);   // Set initial values

    // Operators for the fractional step method
    double* Ru;         // Operator Ru at time n
    double* Rv;         // Operator Rv at time n
    double* Ru_prev;    // Operator Ru at time n-1
    double* Rv_prev;    // Operator Rv at time n-1
    allocateOperatorR(nx, ny, Ru, Rv, Ru_prev, Rv_prev);    // Allocate the previous arrays
    // HEEEEEEEEEEEEEEEEEEEEE
    // computeRu2(Ru_prev, m, u, v, props);                     // Compute initial value of Ru
    // computeRv2(Rv_prev, m, u, v, props);                     // Compute initial value of Rv

    // Predictor velocities
    double* u_pred;     // X-component of predictor velocity
    double* v_pred;     // Y-component of predictor velocity
    allocatePredictorVelocities(nx, ny, u_pred, v_pred);    // Allocate the previous arrays

    // Linear system variables
    double* A;          // Linear system matrix
    double* b;          // Linear system vector
    allocateLinearSystemVariables(nx, ny, A, b);    // Allocate the previous arrays

    std::chrono::steady_clock::time_point begin, end;




    double* mx_staggX = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    if(!mx_staggX) {
        printf("Error: could not allocate enough memory for mx_staggX\n");
        return;
    }

    double* my_staggX = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    if(!my_staggX) {
        printf("Error: could not allocate enough memory for my_staggX\n");
        return;
    }

    double* ux_staggX = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    if(!ux_staggX) {
        printf("Error: could not allocate enough memory for ux_staggX\n");
        return;
    }

    double* uy_staggX = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    if(!uy_staggX) {
        printf("Error: could not allocate enough memory for uy_staggX\n");
        return;
    }


    double* mx_staggY = (double*) calloc((nx+1)*(ny+1), sizeof(double));
    if(!mx_staggY) {
        printf("Error: could not allocate enough memory for mx_staggY\n");
        return;
    }

    double* my_staggY = (double*) calloc((nx+2)*(ny+2), sizeof(double));
    if(!my_staggY) {
        printf("Error: could not allocate enough memory for my_staggY\n");
        return;
    }
    double* vx_staggY = (double*) calloc((nx+1)*(ny+1), sizeof(double));    // ve vw
    if(!vx_staggY) {
        printf("Error: could not allocate enough memory for vx_staggY\n");
        return;
    }
    double* vy_staggY = (double*) calloc((nx+2)*(ny+2), sizeof(double));     // vs vn
    if(!vy_staggY) {
        printf("Error: could not allocate enough memory for vy_staggY\n");
        return;
    }

    // computeMassFlowsStaggX(mx_staggX, my_staggX, m, u, v, props);
    // computeVelocitiesStaggX_CDS(ux_staggX, uy_staggX, nx, ny, u);





    // double* mx_staggY;
    // double* my_staggY;

    // double* Ru_test = (double*) calloc((nx+1)*(ny+2), sizeof(double));
    // if(!Ru_test) {
    //     printf("Error: could not allocate enough memory for Ru_test\n");
    //     return;
    // }
    //
    // // Operator Rv at time n
    // double* Rv_test = (double*) calloc((nx+2)*(ny+1), sizeof(double));
    // if(!Rv_test) {
    //     printf("Error: could not allocate enough memory for Rv_test\n");
    //     return;
    // }

    computeMassFlowsStaggX(mx_staggX, my_staggX, m, u, v, props);
    // computeVelocitiesStaggX_CDS(ux_staggX, uy_staggX, nx, ny, u);
    computeVelocitiesStaggX_QUICK(ux_staggX, uy_staggX, m, mx_staggX, u);
    // computeRu3(Ru_test, m, u, ux_staggX, uy_staggX, mx_staggX, my_staggX, props);
    computeRu3(Ru_prev, m, u, ux_staggX, uy_staggX, mx_staggX, my_staggX, props);

    computeMassFlowsStaggY(mx_staggY, my_staggY, m, u, v, props);
    // computeVelocitiesStaggY_CDS(vx_staggY, vy_staggY, nx, ny, v);
    computeVelocitiesStaggY_QUICK(vx_staggY, vy_staggY, m, my_staggY, v);
    // computeRv3(Rv_test, m, v, vx_staggY, vy_staggY, mx_staggY, my_staggY, props);
    computeRv3(Rv_prev, m, v, vx_staggY, vy_staggY, mx_staggY, my_staggY, props);


    // computeRu2(Ru_prev, m, u, v, props);                     // Compute initial value of Ru
    // computeRv2(Rv_prev, m, u, v, props);                     // Compute initial value of Rv


    // double diffRu = -1;
    // for(int k = 0; k < (nx+1)*(ny+2); k++)
    //     diffRu = std::max(diffRu, std::abs(Ru_prev[k] - Ru_test[k]));
    //
    // double diffRv = -1;
    // for(int k = 0; k < (nx+2)*(ny+1); k++)
    //     diffRv = std::max(diffRv, std::abs(Rv_prev[k] - Rv_test[k]));

    // printf("%20s %10.3e %5s %10.3e\n", "Ru Rv prev diffs", diffRu, "", diffRv);

    // Run the loop until steady state is reached

    double objective_time = 5;

    while(!steady) {

        begin = std::chrono::steady_clock::now();

        // Compute operator R(u) and R(v)
        // computeRu2(Ru_test, m, u, v, props);
        // computeRv2(Rv_test, m, u, v, props);

        computeMassFlowsStaggX(mx_staggX, my_staggX, m, u, v, props);
        // computeVelocitiesStaggX_CDS(ux_staggX, uy_staggX, nx, ny, u);
        computeVelocitiesStaggX_QUICK(ux_staggX, uy_staggX, m, mx_staggX, u);
        computeRu3(Ru, m, u, ux_staggX, uy_staggX, mx_staggX, my_staggX, props);

        computeMassFlowsStaggY(mx_staggY, my_staggY, m, u, v, props);
        // computeVelocitiesStaggY_CDS(vx_staggY, vy_staggY, nx, ny, v);
        computeVelocitiesStaggY_QUICK(vx_staggY, vy_staggY, m, my_staggY, v);
        computeRv3(Rv, m, v, vx_staggY, vy_staggY, mx_staggY, my_staggY, props);

        // // Test
        // computeRu2(Ru_test, m, u, v, props);
        // computeRv2(Rv_test, m, u, v, props);
        //
        // double diffRu = -1;
        // for(int k = 0; k < (nx+1)*(ny+2); k++)
        //     diffRu = std::max(diffRu, std::abs(Ru[k] - Ru_test[k]));
        //
        // double diffRv = -1;
        // for(int k = 0; k < (nx+2)*(ny+1); k++)
        //     diffRv = std::max(diffRv, std::abs(Rv[k] - Rv_test[k]));



        // Compute predictor velocities
        computePredictorVelocityU(u_pred, nx, ny, u, Ru, Ru_prev, props, tstep);
        computePredictorVelocityV(v_pred, nx, ny, v, Rv, Rv_prev, props, tstep);
        lid_driven::setBoundaryPredictorVelocities(u_pred, v_pred, nx, ny, u_ref);
        // Compute discretization coefficients
        computeDiscretizationCoefficients(A, b, m, u_pred, v_pred, props, tstep);
        lid_driven::computeBoundaryDiscretizationCoefficients(A, b, m);
        // Solve linear system
        int exitCodeGS = solveSystemGS(nx+2, ny+2, tol, maxIt, A, b, p);
        if(exitCodeGS == -1) {
            printf("Error solving the linear system. Convergence not achieved.\n");
            return;
        }
        // Compute the velocities at time n+1
        double maxDerivative = 0;   // Max derivative of velocity with respect to time
        computeNextVelocityField(u, v, m, u_pred, v_pred, p, props, tstep, maxDerivative);
        maxDerivative /= tstep;
        // Update current time and number of iterations
        t += tstep;
        it++;
        // Time elapsed
        end = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1e6;


        // printf("%20s %10.3e %5s %10.3e\n", "Ru Rv diffs", diffRu, "", diffRv);
        printf("%6d %5s %10.5f %5s %10.5f %5s %10.5e %5s %10d %5s %.3f\n", it, "", t, "", tstep, "", maxDerivative, "", exitCodeGS, "", elapsed);
        // printf("%20s %10.3e %5s %10.3e\n", "Ru Rv diffs", diffRu, "", diffRv);

        if(t > objective_time) {
            int Re = std::round(props.rho * u_ref * L / props.mu);
            printVariablesToFile(m, u, v, p, Re, 5, std::floor(t));
            objective_time += 5;
        }

        // Check steady state condition
        if(maxDerivative < sstol) {
            steady = true;

            int Re = std::round(props.rho * u_ref * L / props.mu);
            printVariablesToFile(m, u, v, p, Re, 5, std::round(t));

        }
        else {
            // Compute next time step
            computeTimeStep(tstep, m, u, v, props, tol);
            // Update operator R
            std::memcpy(Ru_prev, Ru, (nx+1)*(ny+2)*sizeof(double));
            std::memcpy(Rv_prev, Rv, (nx+2)*(ny+1)*sizeof(double));
        }
    }

    // Free arrays
    free(u);
    free(v);
    free(p);
    free(Ru);
    free(Rv);
    free(Ru_prev);
    free(Rv_prev);
    free(u_pred);
    free(v_pred);
    free(A);
    free(b);

    free(mx_staggX);
    free(my_staggX);
    free(ux_staggX);
    free(uy_staggX);

    free(mx_staggY);
    free(my_staggY);
    free(vx_staggY);
    free(vy_staggY);

    // free(Ru_test);
    // free(Rv_test);
}

void lid_driven::setInitialMaps(double* u, double* v, double* p, const NCMesh m, const double u_ref, const double p_ref) {

    // Mesh size
    int nx = m.getNX(); // Number of control volumes along x axis
    int ny = m.getNY(); // Number of control volumes along y axis

    // Arrays u, v and p allocated using calloc, all bits set to zero by calloc

    // x-velocity (u) initial map
    for(int i = 0; i < nx+1; i++)
        u[(ny+1)*(nx+1)+i] = u_ref;

    // pressure (p) initial map
    for(int k = 0; k < (nx+2)*(ny+2); k++)
        p[k] = p_ref;
}


void lid_driven::setBoundaryPredictorVelocities(double* u_pred, double* v_pred, const int nx, const int ny, const double u_ref) {

    // PREDICTOR VELOCITY X
    // Lower and upper boundaries
    for(int i = 0; i < nx+1; i++) {
        u_pred[i] = 0;              // Lower boundary
        int j = ny + 1;
        u_pred[j*(nx+1)+i] = u_ref; // Upper boundary
    }

    // Left and right boundaries
    for(int j = 1; j < ny+1; j++) {
        u_pred[j*(nx+1)] = 0;           // Left boundary
        u_pred[j*(nx+1)+nx] = 0;        // Right boundary
    }

    // PREDICTOR VELOCITY Y
    // Lower and upper boundaries
    for(int i = 0; i < nx+2; i++) {
        v_pred[i] = 0;                  // Lower boundary
        v_pred[ny*(nx+2)+i] = 0;        // Upper boundary
    }

    // Left and right boundaries
    for(int j = 1; j < ny+1; j++) {
        v_pred[j*(nx+2)] = 0;           // Left boundary
        v_pred[j*(nx+2)+nx+1] = 0;      // Right boundary
    }

}

void lid_driven::computeBoundaryDiscretizationCoefficients(double* A, double* b, const NCMesh m) {

    // Mesh size
    int nx = m.getNX(); // Number of control volumes along x axis
    int ny = m.getNY(); // Number of control volumes along y axis

    // lOWER AND UPPER BOUNDARIES (including corner nodes, i.e. singular nodes)
    for(int i = 0; i < nx+2; i++) {
        // Lower boundary
        int node = i;       // Node number
        A[5*node] = 0;      // South node
        A[5*node+1] = 0;    // West node
        A[5*node+2] = 0;    // East node
        A[5*node+3] = 1;    // North node
        A[5*node+4] = 1;    // Central node
        b[node] = 0;        // Independent term
        // Upper boundary
        node = (ny + 1) * (nx + 2) + i;
        A[5*node] = 1;      // South node
        A[5*node+1] = 0;    // West node
        A[5*node+2] = 0;    // East node
        A[5*node+3] = 0;    // North node
        A[5*node+4] = 1;    // Central node
        b[node] = 0;        // Independent term
    }

    // LEFT AND RIGHT BOUNDARIES
    for(int j = 1; j < ny+1; j++) {
        // Left boundary
        int node = j * (nx + 2);
        A[5*node] = 0;      // South node
        A[5*node+1] = 0;    // West node
        A[5*node+2] = 1;    // East node
        A[5*node+3] = 0;    // North node
        A[5*node+4] = 1;    // Central node
        b[node] = 0;        // Independent term
        // Right boundary
        node = j * (nx + 2) + (nx + 1);
        A[5*node] = 0;      // South node
        A[5*node+1] = 1;    // West node
        A[5*node+2] = 0;    // East node
        A[5*node+3] = 0;    // North node
        A[5*node+4] = 1;    // Central node
        b[node] = 0;        // Independent term
    }
}
