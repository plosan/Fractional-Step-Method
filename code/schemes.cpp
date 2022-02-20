#include "schemes.h"

// Upwind-Difference Scheme
double schemeUDS(double varNeg, double varPos, double m) {
    // Parameters: fix a direction (i, j or k) then
    //  - varNeg    variable var at the node in the negative side of the direction
    //  - varPos    variable var at the node in the positive side of the direction
    //  - m         mass flow through the face
    if(m >= 0)
        return varNeg;
    return varPos;
}

// Central-Difference Scheme
double schemeCDS(double varP, double varF) {
    // Parameters:
    //  - P     Variable var at node P
    //  - C     Variable var at node F
    return 0.5*(varP + varF);
}
