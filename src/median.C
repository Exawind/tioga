#include "median.h"
#include <algorithm>

// Translated verbatim from median.F90
void median(int ix[], double x[], int& n,
                double& xmed) {
  double temp, xhi, xlo, xmax, xmin;
  bool odd;
  int hi, lo, nby2, nby2p1, mid, i, j, k, itemp;

  nby2 = n / 2;
  nby2p1 = nby2 + 1;
  odd = true;

  //     HI & LO are position limits encompassing the median.

  if (n == (2 * nby2)) {odd = false;}
  lo = 1;
  hi = n;
  if (n < 3) {
    if (n < 1) {
      xmed = 0.0;
      return;
    }
    xmed = x[1-1];
    if (n == 1) {return;}
    xmed = 0.5 * (xmed + x[2-1]);
    if (x[2-1] < x[1-1]) {
      temp = x[1-1];
      x[1-1] = x[2-1];
      x[2-1] = temp;
      itemp = ix[1-1];
      ix[1-1] = ix[2-1];
      ix[2-1] = itemp;
    }
    return;
  }

  //     Find median of 1st, middle & last values.
ten:
  mid = (lo + hi) / 2;
  xmed = x[mid-1];
  xlo = x[lo-1];
  xhi = x[hi-1];
  if(xhi < xlo) {
    temp = xhi;
    xhi = xlo;
    xlo = temp;
  }
  if (xmed > xhi) {
    xmed = xhi;
  } else if (xmed < xlo) {
    xmed = xlo;
  }

  // The basic quicksort algorithm to move all values <= the sort key (XMED)
  // to the left-hand end, and all higher values to the other end.
  
  i = lo;
  j = hi;
fifty:
  do {
    if (x[i-1] >= xmed) {break;}
    i = i + 1;
  } while(true);
  do {
    if (x[j-1] <= xmed) {break;}
    j = j - 1;
  } while(true);
  if (i < j) {
    temp = x[i-1];
    x[i-1] = x[j-1];
    x[j-1] = temp;

    itemp=ix[i-1];
    ix[i-1]=ix[j-1];
    ix[j-1]=itemp;

    i = i + 1;
    j = j - 1;

  //     Decide which half the median is in.
    if (i <= j) {goto fifty;}
  }

  if (!odd) {
    if ((j == nby2) && (i == nby2p1)) {goto onehundredthirty;}
    if (j < nby2) {lo = i;}
    if (i > nby2p1) {hi = j;}
    if (i != j) {goto onehundred;}
    if (i == nby2) {lo = nby2;}
    if (j == nby2p1) {hi = nby2p1;}
  } else {
    if (j < nby2p1) {lo = i;}
    if (i > nby2p1) {hi = j;}
    if (i != j) {goto onehundred;}

  // Test whether median has been isolated.

    if (i == nby2p1) {return;}
  }
onehundred:
  if (lo < hi - 1) {goto ten;}

  if (!odd) {
    if (x[nby2p1-1] < x[nby2-1]) {
     temp=x[nby2p1-1];
     x[nby2p1-1]=x[nby2-1];
     x[nby2-1]=temp;
     itemp=ix[nby2p1-1];
     ix[nby2p1-1]=ix[nby2-1];
     ix[nby2-1]=itemp;
    }
    xmed = 0.5*(x[nby2-1] + x[nby2p1-1]);
    return;
  }
  temp = x[lo-1];
  itemp = ix[lo-1];
  if (temp > x[hi-1]) {
    x[lo-1] = x[hi-1];
    x[hi-1] = temp;
    ix[lo-1]=ix[hi-1];
    ix[hi-1]=itemp;
  }
  xmed = x[nby2p1-1];
  return;

  // Special case, N even, J = N/2 & I = J + 1, so the median is
  // between the two halves of the series.   Find max. of the first
  // half & min. of the second half, then average.

onehundredthirty:
  xmax = x[1-1];
  for (k = lo; k <= j; k++) {
    xmax = std::max(xmax, x[k-1]);
  }
  xmin = x[n-1];
  for (k = i; k <= hi; k++) {
    xmin = std::min(xmin, x[k-1]);
  }
  xmed = 0.5*(xmin + xmax);

  return;
}
