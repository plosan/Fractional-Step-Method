#ifndef SCHEMES_H
#define SCHEMES_H

// This header contains the declaration of schemes for the evaluation of convective properties

double schemeUDS(double C, double D, double m); // Upwind Difference Scheme
double schemeCDS(double varP, double varF);     // Centrla Difference Scheme

#endif // SCHEMES_H
