#ifndef SCHEMES_H
#define SCHEMES_H

// This header contains the declaration of schemes for the evaluation of convective properties

// Upwind Difference Scheme
double schemeUDS(double C, double D, double m);

// Central Difference Scheme
double schemeCDS(double varP, double varF);

// Quadratic Upwind Interpolation for Convective Kinematics
double schemeQUICK(double varD, double varU, double varUU, double xD, double xU, double xUU, double xe);

#endif // SCHEMES_H
