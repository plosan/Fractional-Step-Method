#include "solver.h"

int solveSystemGS(const int nx, const int ny, const double tol, const int maxIt, const double* A, const double* b, double* phi) {
    /*
    solveSystemGS: solves the linear system resulting from a 2D convection-diffusion problem in a domain discretized with a cartesian mesh using
    Gauss-Seidel algorithm. It has two criterion to stop the iteration:
        - Let phi* and phi be two consecutive vectors of the sequence produced by Gauss-Seidel algorithm. If the infinity norm of phi-phi* is less
        than tol, then the algorithm stops.
        - If the number of iterations reach maxIt, the algorithm stops.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx        Number of nodes in the X axis                               [const int]
        - ny        Number of nodes in the Y axis                               [const int]
        - tol       Tolerance to stop iteration                                 [const double]
        - maxIt     Maximum number of iterations                                [const int]
        - A         Linear system matrix. Rows: nx*ny, Columns: 5               [double*]
        - b         Vector of indenpendent terms. Rows: nx*ny, Columns: 1       [double*]
        - phi       Solution of the linar system. Rows: nx*ny, Columns: 1       [double*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - phi       Solution of the linar system. Rows: nx*ny, Columns: 1       [double*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Return:
        - int       If the algorithm reaches the max number of iterations with no convergence, it returns -1. Otherwise returns the number of iterations.
    */

    int it = 0;        // Current iteration
    bool convergence = false;   // Boolean variable to tell whether there is convergence or not. False: no convergence, True: convergence
    // Gauss-Seidel iteration
    while(it < maxIt && !convergence) {
        double maxDiff = -1;    // Infinity norm of the difference phi-phi*
        // Lower row nodes
        for(int i = 0; i < nx; i++) {
            int node = i;                                                       // Node whose phi is being computed
            double aux = phi[node];                                             // Previous value of phi[node]
            phi[node] = (b[node] + A[5*node+3] * phi[node+nx]) / A[5*node+4];   // Compute new value
            maxDiff = std::max(maxDiff, std::abs(aux - phi[node]));             // Update infinity norm
        }
        // Central rows nodes
        for(int j = 1; j < ny-1; j++) {
            for(int i = 0; i < nx; i++) {
                int node = j * nx + i;                                          // Node whose phi is being computed
                double aux = phi[node];                                         // Previous value of phi[node]
                phi[node] = (b[node] + A[5*node] * phi[node-nx] + A[5*node+1] * phi[node-1] + A[5*node+2] * phi[node+1] + A[5*node+3] * phi[node+nx]) / A[5*node+4];    // Compute new value
                maxDiff = std::max(maxDiff, std::abs(aux - phi[node]));         // Update infinity norm
            }
        }
        // Upper row nodes
        for(int i = 0; i < nx; i++) {
            int node = (ny - 1) * nx + i;                                       // Node whose phi is being computed
            double aux = phi[node];                                             // Previous value of phi[node]
            phi[node] = (b[node] + A[5*node] * phi[node-nx]) / A[5*node+4];     // Compute new value
            maxDiff = std::max(maxDiff, std::abs(aux - phi[node]));             // Update infinity norm
        }
        // Final checks of the current iteration
        convergence = (maxDiff < tol);  // Convergence condition
        it++;                           // Increase iteration counter
    }
    // Return value
    if(convergence)
        return it;
    return -1;
}
