//
// This file is part of the Tioga software library
//
// Tioga  is a tool for overset grid assembly on parallel distributed systems
// Copyright (C) 2015 Jay Sitaraman
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef TIOGA_UTILS_H
#define TIOGA_UTILS_H

/* header files */
#include "tioga.h"

/* system header files */
#include <iostream>
#include <vector>
#include <unordered_set>

struct Node {
  int id;
  double x,y,z;
  const double eps = 1E-10;
  const double binScaling = 1.0E4; // see note below

  /*< NOTE: Increase this number (e.g.,1.0E5) if lots of x/y/z pts are verry close.
   *
   *  Background Info:
   *    We want to create unique hashes based on coordinates;
   *    by scaling the coordinates, we get better hash spreading.
   *
   *    Example: say we have x-coordinates
   *     x1 = 1.2456
   *     x2 = 1.2467
   *  >> xHash(1E4 * x1) = hash(int(1245.6)) = hash(1245)
   *  >> xHash(1E4 * x2) = hash(int(1246.7)) = hash(1246)
   *     Thus, xHash(x1) != xHash(x2), so no lookup collisions
   *     (so no '==' check).
   *
   *  If we did not perform bin scaling, all coordinates with same whole number
   *  part (left digits of decimal part) would result in the same Hash,
   *  thus requiring the expensive `operator==(const Node& otherNode)' check,
   *  which is an N^2 algorithm (horrific scaling for even a small number of pts).
   */
  Node() { }
  Node(int id,double *geo){
    this->id = id;   /*< node id */
    this->x = geo[0];/*< x-coordinate of point */
    this->y = geo[1];/*< y-coordinate of point */
    this->z = geo[2];/*< z-coordinate of point */
  }

  // checks if two nodes are the same: either by ID or by spatial distance
  bool operator==(const Node& otherNode) const {
    return (otherNode.id == id) || ((abs(this->x - otherNode.x) <= this->eps) &&
                                    (abs(this->y - otherNode.y) <= this->eps) &&
                                    (abs(this->z - otherNode.z) <= this->eps));
  }

  // creates the hash value for inserting unique nodes into an std::unordered_set;
  // if two nodes have the same hash value, the 'operator==' is called.
  struct HashFunction {
    size_t operator()(const Node& node) const {
      size_t xHash = std::hash<int>()(int(node.binScaling*node.x));
      size_t yHash = std::hash<int>()(int(node.binScaling*node.y)) << 1;
      size_t zHash = std::hash<int>()(int(node.binScaling*node.z)) << 2;
      return xHash ^ yHash ^ zHash;
    }
  };
};

/* function declarations */
void findOBB(double *x, double xc[3], double dxc[3], double vec[3][3],
             int nnodes);
int checkHoleMap(double *x, int *nx, int *sam, double *extents);
int checkAdaptiveHoleMap(double *xpt,ADAPTIVE_HOLEMAP *AHM);
void fillHoleMap(int *holeMap, int ix[3], int isym);
void octant_children(uint8_t children_level,uint32_t idx,
                     octant_full_t * q,
                     octant_full_t * c0, octant_full_t * c1,
                     octant_full_t * c2, octant_full_t * c3,
                     octant_full_t * c4, octant_full_t * c5,
                     octant_full_t * c6, octant_full_t * c7);
void octant_children_neighbors(const octant_full_t * q,
                               octant_full_t *c0, octant_full_t *c1,
                               octant_full_t *c2, octant_full_t *c3,
                               octant_full_t *c4, octant_full_t *c5,
                               octant_full_t *c6, octant_full_t *c7);
void floodfill_level(level_octant_t *level);
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

void qcoord_to_vertex(qcoord_t x,qcoord_t y,qcoord_t z,double *vertices,double vxyz[3]);
char checkFaceBoundaryNodes(int *nodes,const char *bcnodeflag,const int numfaceverts,const int *faceConn,const char *duplicatenodeflag);
int triBoxOverlap(double boxcenter[3],double boxhalfsize[3],double *pt1,double *pt2,double *pt3);

/* inline functions */
//#include <malloc.h>
//static inline
//double memory_usage(int mpi_rank,int timestep,int display){
//
//    /* get malloc info structure */
//    struct mallinfo my_mallinfo = mallinfo();
//
//    /*total memory reserved by the system for malloc currently */
//    double reserved_mem = my_mallinfo.arena;
//
//    /* get all the memory currently allocated to user by malloc, etc. */
//    double used_mem = my_mallinfo.hblkhd
//                    + my_mallinfo.usmblks
//                    + my_mallinfo.uordblks;
//
//    /* get memory not currently allocated to user but malloc controls */
//    double free_mem = my_mallinfo.fsmblks
//                    + my_mallinfo.fordblks;
//
//    /* get number of items currently allocated */
//    /* double number_allocated = my_mallinfo.ordblks + my_mallinfo.smblks; */
//
//    /* Print out concise malloc info line */
//    if(display && mpi_rank == 0){
//        printf("Step[%d]: %f MB(%.0f) malloc: %f MB reserved (%.0f unused)\n",
//            timestep,
//            used_mem / (1024.0 * 1024.0),
//            used_mem,
//            reserved_mem / (1024.0 * 1024.0),
//            free_mem);
//
//        if(mpi_rank == 0){
//            FILE *fp;
//            char filename[] = "tiogaMemUsage.dat";
//            fp=fopen(filename,"a");
//            fprintf(fp,"Step[%d]: %f MB(%.0f) malloc: %f MB reserved (%.0f unused)\n",
//                    timestep,used_mem / (1024.0 * 1024.0),used_mem,reserved_mem / (1024.0 * 1024.0), free_mem);
//            fclose(fp);
//        }
//
//    }
//    return used_mem;
//}

#endif /* TIOGA_UTILS_H */
