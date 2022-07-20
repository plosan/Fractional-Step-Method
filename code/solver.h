#ifndef SOLVER_H
#define SOLVER_H

#include <algorithm>

int solveSystemGS(const int nx, const int ny, const double tol, const int maxIt, const double* A, const double* b, double* phi);

#endif // SOLVER_H
