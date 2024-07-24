#ifndef TIOGA_MATH_H
#define TIOGA_MATH_H
double computeCellVolume(double xv[8][3], int nvert);
double tdot_product(double a[3], double b[3], double c[3]);
void computeNodalWeights(
    double xv[8][3], double* xp, double frac[8], int nvert);
#endif
