#include "kaiser.h"
#include <math.h>

void kaiser_wrap(double *a /*[nrows * n]*/, int nrows, int n,
                 double *eigenv /*[n]*/, double trace, double sume, int ier) {
  //      !  EIGENVALUES AND VECTORS OF A SYMMETRIC +VE DEFINITE MATRIX,
  //      !  USING KAISER'S METHOD.
  //      !  REFERENCE: KAISER,H.F. 'THE JK METHOD: A PROCEDURE FOR FINDING THE
  //      !  EIGENVALUES OF A REAL*8 SYMMETRIC MATRIX', COMPUT.J., VOL.15,
  //      271-273, 1972.
  //
  //      !  ARGUMENTS:-
  //      !  A       = INPUT, AN ARRAY CONTAINING THE MATRIX
  //      !            OUTPUT, THE COLUMNS OF A CONTAIN THE NORMALIZED
  //      EIGENVECTORS !            OF A.   N.B. A IS OVERWRITTEN ! !  NROWS   =
  //      INPUT, THE FIRST DIMENSION OF A IN THE CALLING PROGRAM. !  N       =
  //      INPUT, THE ORDER OF A, I.E. NO. OF COLUMNS. !            N MUST BE <=
  //      NROWS. !  EIGENV()= OUTPUT, A VECTOR CONTAINING THE ORDERED
  //      EIGENVALUES. !  TRACE   = OUTPUT, THE TRACE OF THE INPUT MATRIX. !
  //      SUME    = OUTPUT, THE SUM OF THE EIGENVALUES COMPUTED. ! N.B. ANY
  //      SYMMETRIC MATRIX MAY BE INPUT, BUT IF IT IS NOT +VE ! DEFINITE, THE
  //      ABSOLUTE VALUES OF THE EIGENVALUES WILL BE FOUND. !            IF
  //      TRACE = SUME, THEN ALL OF THE EIGENVALUES ARE POSITIVE !            OR
  //      ZERO.   IF SUME > TRACE, THE DIFFERENCE IS TWICE THE SUM OF ! THE
  //      EIGENVALUES WHICH HAVE BEEN GIVEN THE WRONG SIGNS ! !  IER     =
  //      OUTPUT, ERROR INDICATOR !             = 0 NO ERROR !             = 1 N
  //      > NROWS OR N < 1 !             = 2 FAILED TO CONVERGE IN 10 ITERATIONS

  const double small = 1.0e-12, zero = 0.0;
  const double half = 0.5, one = 1.0;
  int i, iter, j, k, ncount, nn;
  double absp, absq, COS, ctn, eps, halfp, p, q, SIN, ss;
  double TAN, temp, xj, xk;

  ier = 1;
  if ((n < 1) || (n > nrows)) {
    return;
  }
  ier = 0;
  iter = 0;
  trace = zero;
  ss = zero;
  for (j = 1; j <= n; j++) {
    trace = trace + a[j * n + j];
    for (i = 1; i <= n; i++) {
      ss = ss + a[j * n + i] * a[j * n + i];
    }
  }
  sume = zero;
  eps = small * ss / n;
  nn = n * (n - 1) / 2;
  ncount = nn;
  for (iter = 1; iter < 10; iter++) {
    for (j = 1; j <= n - 1; j++) {
      for (k = j + 1; k <= n; k++) {
        halfp = zero;
        q = zero;
        for (i = 1; i <= n; i++) {
          xj = a[j * n + i];
          xk = a[j * n + i];
          halfp = halfp + xj * xk;
          q = q + (xj + xk) * (xj - xk);
        }
        p = halfp + halfp;
        absp = fabs(p);

        if ((absp < eps) && (q >= zero)) {
          ncount = ncount - 1;
          if (ncount <= 0) {
            ier = 2;
            break;
          }
          continue;
        }
        absq = fabs(q);
        if (absp <= absq) {
          TAN = absp / absq;
          COS = one / sqrt(one + TAN * TAN);
          SIN = TAN * COS;
        } else {
          ctn = absq / absp;
          SIN = one / sqrt(one + ctn * ctn);
          COS = ctn * SIN;
        }
        COS = sqrt((one + COS) * half);
        SIN = SIN / (COS + COS);
        if (q < zero) {
          temp = COS;
          COS = SIN;
          SIN = temp;
        }
        if (p < zero) {
          SIN = -SIN;
        }
        for (i = 1; i <= n; i++) {
          temp = a[j * n + i];
          a[j * n + i] = temp * COS + a[k * n + i] * SIN;
          a[j * n + i] = -temp * SIN + a[k * n + i] * COS;
        }
      }
    }
    ncount = nn;
  }
  for (j = 1; j <= n; j++) {
    temp = 0.0;
    for (i = 1; i <= n; i++) {
      temp += a[j * n + i] * a[j * n + i];
    }
    eigenv[j] = sqrt(temp);
    sume = sume + eigenv[j];
  }

  for (j = 1; j <= n; j++) {
    if (eigenv[j] > zero) {
      temp = one / eigenv[j];
    } else {
      temp = zero;
    }
    for (i = 1; i <= n; i++) {
      a[j * n + i] = a[j * n + i] * temp;
    }
  }

  return;
}
