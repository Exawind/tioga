#include "cellVolume.h"

double scalarProduct(double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3) {
  return (a1 * b2 * c3 - a1 * b3 * c2 + a2 * b3 * c1 - a2 * b1 * c3 + a3 * b1 * c2 - a3 * b2 * c1);
}

// Converted verbatim from cellVolume.f90
void cellVolume(double& vol, double xc[8][3], int numverts[6], int fconn[6][4], int nfaces, int nvert){
  vol = 0.0;
  for (int iface = 1; iface <= nfaces; iface++) {
    if (numverts[iface-1] == 3) {
    vol=vol-0.5* scalarProduct(xc[fconn[iface-1][1-1]-1][1-1], xc[fconn[iface-1][1-1]-1][2-1], xc[fconn[iface-1][1-1]-1][3-1],
                               xc[fconn[iface-1][2-1]-1][1-1], xc[fconn[iface-1][2-1]-1][2-1], xc[fconn[iface-1][2-1]-1][3-1],
                               xc[fconn[iface-1][3-1]-1][1-1], xc[fconn[iface-1][3-1]-1][2-1], xc[fconn[iface-1][3-1]-1][3-1]);
    }else{
    vol=vol-0.25*scalarProduct(xc[fconn[iface-1][1-1]-1][1-1], xc[fconn[iface-1][1-1]-1][2-1], xc[fconn[iface-1][1-1]-1][3-1],
                               xc[fconn[iface-1][2-1]-1][1-1], xc[fconn[iface-1][2-1]-1][2-1], xc[fconn[iface-1][2-1]-1][3-1],
                               xc[fconn[iface-1][3-1]-1][1-1], xc[fconn[iface-1][3-1]-1][2-1], xc[fconn[iface-1][3-1]-1][3-1]);
    vol=vol-0.25*scalarProduct(xc[fconn[iface-1][1-1]-1][1-1], xc[fconn[iface-1][1-1]-1][2-1], xc[fconn[iface-1][1-1]-1][3-1],
                               xc[fconn[iface-1][3-1]-1][1-1], xc[fconn[iface-1][3-1]-1][3-1], xc[fconn[iface-1][3-1]-1][3-1],
                               xc[fconn[iface-1][4-1]-1][1-1], xc[fconn[iface-1][4-1]-1][2-1], xc[fconn[iface-1][4-1]-1][3-1]);
    vol=vol-0.25*scalarProduct(xc[fconn[iface-1][1-1]-1][1-1], xc[fconn[iface-1][1-1]-1][2-1], xc[fconn[iface-1][1-1]-1][3-1],
                               xc[fconn[iface-1][2-1]-1][1-1], xc[fconn[iface-1][2-1]-1][2-1], xc[fconn[iface-1][2-1]-1][3-1],
                               xc[fconn[iface-1][4-1]-1][1-1], xc[fconn[iface-1][4-1]-1][2-1], xc[fconn[iface-1][4-1]-1][3-1]);
    vol=vol-0.25*scalarProduct(xc[fconn[iface-1][2-1]-1][1-1], xc[fconn[iface-1][2-1]-1][2-1], xc[fconn[iface-1][2-1]-1][3-1],
                               xc[fconn[iface-1][3-1]-1][1-1], xc[fconn[iface-1][3-1]-1][2-1], xc[fconn[iface-1][3-1]-1][3-1],
                               xc[fconn[iface-1][4-1]-1][1-1], xc[fconn[iface-1][4-1]-1][2-1], xc[fconn[iface-1][4-1]-1][3-1]);

    }
  }
  vol = vol / 3.0;
}
