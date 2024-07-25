#ifndef TIOGA_MATH_H
#define TIOGA_MATH_H
double computeCellVolume(double xv[8][3], int nvert);
double tdot_product(const double a[3], const double b[3], const double c[3]);
void computeNodalWeights(
    double xv[8][3], const double* xp, double frac[8], int nvert);
#endif
