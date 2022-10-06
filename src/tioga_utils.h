#ifndef TIOGA_UTILS_H
#define TIOGA_UTILS_H
void findOBB(double *x, double xc[3], double dxc[3], double vec[3][3],
             int nnodes);
int checkHoleMap(double *x, int *nx, int *sam, double *extents);
void fillHoleMap(int *holeMap, int ix[3], int isym);
int obbIntersectCheck(double vA[3][3], double xA[3], double dxA[3],
                      double vB[3][3], double xB[3], double dxB[3]);
void getobbcoords(double xc[3], double dxc[3], double vec[3][3],
                  double xv[8][3]);
void transform2OBB(double xv[3], double xc[3], double vec[3][3], double xd[3]);
void writebbox(OBB *obb, int bid);
void writebboxdiv(OBB *obb, int bid);
void writePoints(double *x, int nsearch, int bid);
void uniquenodes(double *x, int *meshtag, double *rtag, int *itag, int *nn);
void uniqNodesTree(double *coord, int *itag, double *rtag, int *meshtag,
                   int *elementsAvailable, int ndim, int nav);
void uniquenodes_octree(double *x, int *meshtag, double *rtag, int *itag,
                        int *nn);
#endif
