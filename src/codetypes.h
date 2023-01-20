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

#ifndef CODETYPES_H
#define CODETYPES_H

//#define MPICH_SKIP_MPICXX
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include "mpi.h"
/*====================================================================*/
/*  Floating point definition                                         */
/*====================================================================*/
# define REAL double
typedef int32_t qcoord_t;
/*====================================================================*/
/*  Base for indexing (0 or 1) */
/*====================================================================*/
# define BASE 1

/*====================================================================*/
/*  Define arithmetic constants                                       */
/*====================================================================*/
// #define ZERO               0.0e+00
// #define ONE                1.0e+00
// #define TWO                2.0e+00
// #define THREE              3.0e+00
// #define FOUR               4.0e+00
// #define HALF               0.5e+00
// #define THIRD              0.333333333e+00
// #define FIFTH              0.2
// #define PI                 3.1415926535897932e+00
// #define RAD2DEG            (180.0/PI)
// #define DEG2RAD            (PI/180.0)
#define BIGVALUE           1.0e+15
#define BIGINT             2147483647
#define TOL                1.0e-10
#define HOLEMAPSIZE        192
// #define NFRINGE            3
// #define NVAR               6
#define WALLNODETYPE       0
#define OUTERNODETYPE      1
/*==================================================================*/
/* ADAPTIVE HOLE MAP OCTANT INFO                                    */
/*==================================================================*/
#define USE_ADAPTIVE_HOLEMAP 1 // [0] original hole map, [1] adaptive hole map
#define INTERSECT_ALG        1 // [0] point-box inclusion only
                               // [1] face-box intersection (water-tight)

/* Fixed Octree Constraints: Do Not Change */
#define OCTANT_MAXLEVEL     19
#define OCTANT_CHILDREN     8

/** The length of a side of the root quadrant */
#define OCTANT_ROOT_LEN   ((qcoord_t) 1 << OCTANT_MAXLEVEL)
/** The length of a octant of level l */
#define OCTANT_LEN(l) ((qcoord_t) 1 << (OCTANT_MAXLEVEL - (l)))
/** Conversion from integer coordinates to double coordinates */
#define INT2DBL ((double) 1.0 / (double) OCTANT_ROOT_LEN)

#define OUTSIDE_SB 0
#define INSIDE_SB  1
#define WALL_SB    2
/*==================================================================*/
/* inline debugging tools                                           */
/*==================================================================*/
# define TRACEI(x)  printf("#tioga:\t"#x" =%d\n",x);
# define TRACED(x)  printf("#tioga:\t"#x" =%.16e\n",x);
# define TIOGA_MIN(x,y)  (x) < (y) ? (x) : (y)
# define TIOGA_MAX(x,y)  (x) > (y) ? (x) : (y)
# define TIOGA_FREE(a1)  {free(a1);a1=NULL;}
// # define debug(x,y)  printf("#tioga:\t"#x"=%d,"#y"=%d\n",x,y);
// # define stdwrite(x) if (myid==0) printf("#tioga:\t"#x"\n");
// # define dstr(x) printf("#tioga:\t"#x"\n");
// # define ditch(x,y) {dstr(x); TRACEI(y); MPI_Abort(MPI_COMM_WORLD,ierr);}
/*====================================================================*/
/*  Numerical Tools                                                   */
/*====================================================================*/
// #define Sign(a1,a2) (((a2) < ZERO)? - fabs(a1): fabs(a1))
#define TIOGA_Max(a1,a2) (((a1) >= (a2))? (a1): (a2))
#define TIOGA_Min(a1,a2) (((a1) <= (a2))? (a1): (a2))
// #define Abs(aa) (((aa) >= 0)? (aa): -(aa))
// #define Round(aa) (int) ((fabs((aa) - floor(aa)) >= HALF)? ceil(aa): floor(aa))
// #define swap(a,b) { a=a+b;b=a-b;a=a-b;}
/*===================================================================*/
/* Code specific types                                               */
/*===================================================================*/
#define XLO 0
#define XHI 1
#define YLO 2
#define YHI 3
#define ZLO 4
#define ZHI 5

typedef struct {
  double lo;  /**< lower bound */
  double hi;  /**< upper bound */
} bound_t;

typedef struct {
  bound_t x;  /**< x bounds */
  bound_t y;  /**< y bounds */
  bound_t z;  /**< z bounds */
} box_t;

/** The 3D full octant datatype: 130 bytes per octant */
typedef struct octant_full
{
  qcoord_t x, y, z; /**< [12B] binary coordinates */
  uint32_t id;      /**< [4B] element id on level */
  uint8_t filltype; /**< [1B] floodfill: [0] inside SB, [1] outside SB, [2] hole SB */
  uint8_t refined;  /**< [1B] flag if refined (i.e. is a parent) */
  struct octant_full *nhbr[6];     /**< [48B] neighbor octant list */
  struct octant_full *children[8]; /**< [64B] children octants if refined */
} octant_full_t;

/** 3D octant datatype: 48 bytes per octant */
typedef struct octant
{
  qcoord_t x, y, z;     /**< [12B] binary coordinates */
  uint8_t filltype;     /**< [1B]  floodfill: [0] inside SB, [1] outside SB, [2] hole SB */
  uint8_t leafflag;     /**< [1B]  flag if refined (i.e. is a parent) */
  uint8_t pad[2];       /**< [2B]  padding */
  uint32_t children[8]; /**< [32B] children octant IDs */
} octant_t;

typedef struct level_octant
{
  uint32_t elem_count;      /**< number of octants in level */
  uint8_t level_id;         /**< level number */
  std::vector<octant_full_t> octants; /**< [elem_count] locally stored octants */
} level_octant_t;

typedef struct level
{
  uint8_t level_id;     /**< level number */
  uint32_t elem_count;  /**< number of octants in level */
  std::vector<octant_t> octants; /**< [elem_count] octant list */
} level_t;

typedef struct ADAPTIVE_HOLEMAP_OCTANT
{
  int8_t existWall;     /**< flag to indicate map contains wall */
  double extents_lo[3]; /**< lower coordinates of tree */
  double extents_hi[3]; /**< upper coordinates of tree */

  uint8_t nlevel;       /**< number of levels */
  level_octant_t levels[OCTANT_MAXLEVEL];
} ADAPTIVE_HOLEMAP_OCTANT;

typedef struct {
  uint8_t nlevel;       /**< number of levels in map */
  double extents_lo[3]; /**< lower coordinates of tree */
  double extents_hi[3]; /**< upper coordinates of tree */
  uint64_t leaf_count;  /**< total leaf octant count */
  uint64_t elem_count;  /**< total octant count */
} ahm_meta_t ;

typedef struct ADAPTIVE_HOLEMAP
{
  uint8_t existWall;    /**< flag to indicate map contains wall */
  ahm_meta_t meta;      /**< adaptive hole map meta data */
  level_t levels[OCTANT_MAXLEVEL];
} ADAPTIVE_HOLEMAP;

typedef struct HOLEMAP
{
  int existWall;
  int nx[3];
  int *samLocal;
  int *sam;
  double extents[6];
} HOLEMAP;

typedef struct OBB
{
  double xc[3];
  double dxc[3];
  double vec[3][3];

  int comm_idx;    /* Index in comm map for this OBB                       */
  int iblk_local;  /* Index of this mesh block                             */
  int iblk_remote; /* Index of the remote mesh block (intersecting  pair)  */
  int tag_remote;
  int send_tag;
  int recv_tag;
} OBB;

typedef struct DONORLIST
{
  int donorData[4];
  double donorRes;
  double receptorRes;
  int cancel;
  struct DONORLIST *next;
} DONORLIST;

typedef struct PACKET
{
  int nints;
  int nreals;
  int *intData;
  REAL *realData;
} PACKET;

typedef struct INTERPLIST
{
  int cancel;
  int nweights;
  int receptorInfo[3];
  double xtmp[3];
  int *inode;
  double *weights;
} INTERPLIST;

typedef struct INTERPLIST2
{
  int cancel;
  int nweights;
  int receptorInfo[3];
  int *inode;
  double *weights;
  struct INTERPLIST2 *next;
} INTERPLIST2;

typedef struct INTEGERLIST
{
  int inode;
  struct INTEGERLIST *next;
} INTEGERLIST;

typedef struct INTEGERLIST2
{
  int intDataSize,realDataSize;
  int *intData;
  double *realData;
  struct INTEGERLIST2 *next;
} INTEGERLIST2;

#endif /* CODETYPES_H */