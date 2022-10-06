#include "kaiser.h"
#include <math.h>

void kaiser_cpp(double a[3][3], int nrows, int n, double eigenv[3],
                double trace, double sume) {
  const double small = 1.0e-12;
  const double zero = 0.0;
  const double half = 0.5;
  const double one = 1.0;
  int i, iter, j, k, ncount, nn;
  double absp, absq, COS, ctn, eps, halfp, p, q, SIN, ss;
  double TAN, temp, xj, xk;
  int ier;

  ier = 1;
  if ((n < 1) || (n > nrows)) {
    return;
  }
  ier = 0;
  iter = 0;
  trace = zero;
  ss = zero;
  for (j = 1; j <= n; j++) {
    trace = trace + a[j - 1][j - 1];
    for (i = 1; i <= n; i++) {
      ss = ss + a[j - 1][i - 1] * a[j - 1][i - 1];
    }
  }
  sume = zero;
  eps = small * ss / n;
  nn = n * (n - 1) / 2;
  ncount = nn;

  //   ORTHOGONALIZE PAIRS OF COLUMNS J & K, K > J.
twenty:
  for (j = 1; j <= (n - 1); j++) {
    for (k = j + 1; k <= n; k++) {

      //   CALCULATE PLANAR ROTATION REQUIRED

      halfp = zero;
      q = zero;
      for (i = 1; i <= n; i++) {
        xj = a[j - 1][i - 1];
        xk = a[k - 1][i - 1];
        halfp = halfp + xj * xk;
        q = q + (xj + xk) * (xj - xk);
      }
      p = halfp + halfp;
      absp = std::fabs(p);

      //   If P is very small, the vectors are almost orthogonal.
      //   Skip the rotation if Q >= 0 (correct ordering).

      if ((absp < eps) && (q >= zero)) {
        ncount = ncount - 1;
        if (ncount <= 0) {
          goto onehundredsixty;
        }
        continue;
      }

      //   Rotation needed.

      absq = std::fabs(q);
      if (absp <= absq) {
        TAN = absp / absq;
        COS = one / std::sqrt(one + TAN * TAN);
        SIN = TAN * COS;
      } else {
        ctn = absq / absp;
        SIN = one / std::sqrt(one + ctn * ctn);
        COS = ctn * SIN;
      }
      COS = std::sqrt((one + COS) * half);
      SIN = SIN / (COS + COS);
      if (q < zero) {
        temp = COS;
        COS = SIN;
        SIN = temp;
      }
      if (p < zero) {
        SIN = -SIN;
      }

      //   PERFORM ROTATION

      for (i = 1; i <= n; i++) {
        temp = a[j - 1][i - 1];
        a[j - 1][i - 1] = temp * COS + a[k - 1][i - 1] * SIN;
        a[k - 1][i - 1] = -temp * SIN + a[k - 1][i - 1] * COS;
      }
    }
  }
  ncount = nn;
  iter = iter + 1;
  if (iter < 10) {
    goto twenty;
  }
  ier = 2;

  //   CONVERGED, OR GAVE UP AFTER 10 ITERATIONS
onehundredsixty:
  for (j = 1; j <= n; j++) {
    temp = 0.0;
    for (int m = 1; m <= n; m++) {
      temp += (a[j - 1][m - 1] * a[j - 1][m - 1]);
    }
    eigenv[j - 1] = std::sqrt(temp);
    sume = sume + eigenv[j - 1];
  }

  //   SCALE COLUMNS TO HAVE UNIT LENGTH

  for (j = 1; j <= n; j++) {
    if (eigenv[j - 1] > zero) {
      temp = one / eigenv[j - 1];
    } else {
      temp = zero;
    }
    for (int m = 1; m <= n; m++) {
      a[j - 1][m - 1] = a[j - 1][m - 1] * temp;
    }
  }

  return;
}
