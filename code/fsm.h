// fsm.h
#ifndef FSM_H
#define FSM_H

#include "NCMesh.h"

// Data types
struct Properties {
    double rho;
    double mu;
};

// Variable allocation
void allocateFluidVariables(const int nx, const int ny, double* &u, double* &v, double* &p);
void allocateOperatorR(const int nx, const int ny, double* &Ru, double* &Rv, double* &Ru_prev, double* &Rv_prev);
void allocatePredictorVelocities(const int nx, const int ny, double* &u_pred, double* &v_pred);
void allocateLinearSystemVariables(const int nx, const int ny, double* &A, double* &b);

// Operator Ru and Rv
void computeMassFlowsStaggX(double* mx, double* my, const NCMesh m, const double* u, const double* v);
void computeMassFlowsStaggY(double* mx, double* my, const NCMesh m, const double* u, const double* v);

void computeRu(double* Ru, const NCMesh m, const double* u, const double* v, const Properties props);
void computeRv(double* Rv, const NCMesh m, const double* u, const double* v, const Properties props);

void computePredictorVelocityU(double* u_pred, const int nx, const int ny, const double* u, const double* Ru, const double* Ru_prev, const Properties props, const double tstep);
void computePredictorVelocityV(double* v_pred, const int nx, const int ny, const double* v, const double* Rv, const double* Rv_prev, const Properties props, const double tstep);

void computeDiscretizationCoefficients(double* A, double* b, const NCMesh m, const double* u_pred, const double* v_pred, const Properties props, const double tstep);

void computeVelocity(double* u, double* v, const NCMesh m, const double* u_pred, const double* v_pred, const double* p, const Properties props, const double tstep, double& maxDiff);
void computeVelocityU(double* u, const NCMesh m, const double* u_pred, const double* p, const Properties props, const double tstep, double& maxDiff);
void computeVelocityV(double* v, const NCMesh m, const double* v_pred, const double* p, const Properties props, const double tstep, double& maxDiff);
void computeTimeStep(double &tstep, const NCMesh m, const double* u, const double* v, const Properties props, const double zeroTol);
void updateOperatorR(double* Ru_prev, double* Rv_prev, const double* Ru, const double* Rv, const int nx, const int ny);

void computeCenteredNodesVelocities(double* u_col, double* v_col, const int nx, const int ny, const double* u, const double* v);
void printVelocityToFile(const NCMesh m, double* u_col, double* v_col, const char* filename, const int precision);
void printVelocityUToFile(const NCMesh m, double* u_col, const char* filename, const int precision);
void printVelocityVToFile(const NCMesh m, double* v_col, const char* filename, const int precision);
void printPressureToFile(const NCMesh m, double* p, const char* filename, const int precision);

// Lid-driven cavity
namespace lid_driven {
    void mainLoop(const double rho, const double mu, const double u_ref, const double p_ref, const double L, const int nx, const int ny, const double tstep, const int maxIt, const double tol, const double sstol);
    void setInitialMaps(double* u, double* v, double* p, const NCMesh m, const double u_ref, const double p_ref);
    void setBoundaryPredictorVelocities(double* u_pred, double* v_pred, const int nx, const int ny, const double u_ref);
    void computeBoundaryDiscretizationCoefficients(double* A, double* b, const NCMesh m);
};

#endif // FSM_H
