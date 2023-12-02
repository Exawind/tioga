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

/* system header files */
#include <limits>
#include <vector>
#include <array>
#include <cstring>

/* header files */
#include "codetypes.h"
#include "tioga.h"
#include "tioga_utils.h"

using namespace TIOGA;

/**
 * Create hole maps for all grids
 * this routine is not efficient
 * since it does mutiple global reduce ops
 * have to change it at a later date when
 * there is more time to develop code
 */
void tioga::getHoleMap(void)
{
  int i,j,k,m;
  int ii,jj,kk;
  // double wbox[6];
  std::vector<std::array<double,6>> wbox(nblocks);
  //int existWall;
  std::vector<int> existWall(nblocks);
  int meshtag,maxtag, mtagtmp;
  int *existHoleLocal;
  int *existHole;
  double *bboxLocal;
  double *bboxGlobal;
  double ds[3],dsmax,dsbox;
  int bufferSize;
  FILE *fp;
  char fname[80];
  char intstring[12];
 //
 // get the local bounding box
 //
  meshtag = -BIGINT; //std::numeric_limits<int>::lowest();
  for (int i=0; i<nblocks; i++) {
    auto& mb = mblocks[i];
    mb->getWallBounds(&mtagtmp,&existWall[i],wbox[i].data());
    if (mtagtmp > meshtag) meshtag = mtagtmp;
  }
  MPI_Allreduce(&meshtag,&maxtag,1,MPI_INT,MPI_MAX,scomm);
 //
 if (holeMap)
   {
     for(i=0;i<nmesh;i++)
       if (holeMap[i].existWall) TIOGA_FREE(holeMap[i].sam);
     delete [] holeMap;
   }
 holeMap=new HOLEMAP[maxtag];
 //
 existHoleLocal=(int *)malloc(sizeof(int)*maxtag);
 existHole=(int *)malloc(sizeof(int)*maxtag);
 //
 for(i=0;i<maxtag;i++) existHole[i]=existHoleLocal[i]=0;
 //
 for (int i=0; i<nblocks; i++) {
   existHoleLocal[mtags[i]-1]=existWall[i];
 }
 //
 MPI_Allreduce(existHoleLocal,existHole,maxtag,MPI_INT,MPI_MAX,scomm);
 //
 for(i=0;i<maxtag;i++) holeMap[i].existWall=existHole[i];
 //
 bboxLocal=(double *) malloc(sizeof(double)*6*maxtag);
 bboxGlobal=(double *) malloc(sizeof(double)*6*maxtag);
 //
 for(i=0;i<3*maxtag;i++) bboxLocal[i]=BIGVALUE;
 for(i=0;i<3*maxtag;i++) bboxLocal[i+3*maxtag]=-BIGVALUE;
 for(i=0;i<3*maxtag;i++) bboxGlobal[i]=BIGVALUE;
 for(i=0;i<3*maxtag;i++) bboxGlobal[i+3*maxtag]=-BIGVALUE;

 //
 for (int n=0; n<nblocks; n++) {
   meshtag = mtags[n];
   for (i = 0; i < 3; i++) {
     bboxLocal[3 * (meshtag - 1) + i] = wbox[n][i];
     bboxLocal[3 * (meshtag - 1) + i + 3 * maxtag] = wbox[n][i + 3];
   }
 }
 //
 // get the global bounding box info across all the
 // partitions for all meshes
 //
 MPI_Allreduce(bboxLocal, bboxGlobal, 3 * maxtag, MPI_DOUBLE, MPI_MIN, scomm);
 MPI_Allreduce(
   &(bboxLocal[3 * maxtag]), &(bboxGlobal[3 * maxtag]), 3 * maxtag, MPI_DOUBLE,
   MPI_MAX, scomm);
 //
 // find the bounding box for each mesh
 // from the globally reduced data
 //
 for (i = 0; i < maxtag; i++) {
   if (holeMap[i].existWall) {
     for (j = 0; j < 3; j++) {
       holeMap[i].extents[j] = bboxGlobal[3 * i + j];
       holeMap[i].extents[j + 3] = bboxGlobal[3 * i + j + 3 * maxtag];
       ds[j] = holeMap[i].extents[j + 3] - holeMap[i].extents[j];
     }
     dsmax = std::max(ds[0], ds[1]);
     dsmax = std::max(dsmax, ds[2]);
     dsbox = dsmax / HOLEMAPSIZE;

     for (j = 0; j < 3; j++) {
       holeMap[i].extents[j] -= (2 * dsbox);
       holeMap[i].extents[j + 3] += (2 * dsbox);
       holeMap[i].nx[j] = floor(
         std::max((holeMap[i].extents[j + 3] - holeMap[i].extents[j]) / dsbox, 1.0));
     }
     bufferSize = holeMap[i].nx[0] * holeMap[i].nx[1] * holeMap[i].nx[2];
     holeMap[i].sam = (int*)malloc(sizeof(int) * bufferSize);
     holeMap[i].samLocal = (int*)malloc(sizeof(int) * bufferSize);
     for (j = 0; j < bufferSize; j++)
       holeMap[i].sam[j] = holeMap[i].samLocal[j] = 0;
   }
 }
 //
 // mark the wall boundary cells in the holeMap
 //
 for (int ib=0;ib<nblocks;ib++) {
   auto& mb = mblocks[ib];
   meshtag = mb->getMeshTag();
   if (holeMap[meshtag - 1].existWall) {
    mb->markWallBoundary(
      holeMap[meshtag - 1].samLocal, holeMap[meshtag - 1].nx,
      holeMap[meshtag - 1].extents);
   }
 }
 //
 // allreduce the holeMap of each mesh
 //
 for(i=0;i<maxtag;i++)
   {
    if (holeMap[i].existWall)
     {
      bufferSize=holeMap[i].nx[0]*holeMap[i].nx[1]*holeMap[i].nx[2];
      MPI_Allreduce(holeMap[i].samLocal,holeMap[i].sam,bufferSize,MPI_INT,MPI_MAX,scomm);
     }
   }
 //
 for(i=0;i<maxtag;i++)
   if (holeMap[i].existWall) TIOGA_FREE(holeMap[i].samLocal);
 //
 // set the global number of meshes to maxtag
 //
 nmesh=maxtag;
 //
 // now fill the holeMap
 //
 for(i=0;i<maxtag;i++)
   if (holeMap[i].existWall) fillHoleMap(holeMap[i].sam,holeMap[i].nx,isym);
 //
 // output the hole map
 //
 //this->outputHoleMap();
 //
 // free local memory
 //
 TIOGA_FREE(existHoleLocal);
 TIOGA_FREE(existHole);
 TIOGA_FREE(bboxLocal);
 TIOGA_FREE(bboxGlobal);
}

/**
 * Create adaptive hole maps for all grids
 * this routine is not efficient
 * since it does multiple global reduce ops
 * have to change it at a later date when
 * there is more time to develop code
 */
void tioga::getAdaptiveHoleMap(void){
  int maxtag,maxtagLocal;
  int mbi,mi;
  int level_id;

  /* =========================== */
  /* A: count max number of tags */
  /* =========================== */
  maxtagLocal = -BIGINT;
  for(mbi=0; mbi<nblocks; mbi++){
    auto& mb = mblocks[mbi];
    int mbtag = mb->getMeshTag();
    maxtagLocal = (maxtagLocal < mbtag) ? mbtag:maxtagLocal;
  }
  MPI_Allreduce(&maxtagLocal,&maxtag,1,MPI_INT,MPI_MAX,scomm);

  /* =============================== */
  /* B: reallocate adaptive hole map */
  /* =============================== */
  if(adaptiveHoleMap) delete[] adaptiveHoleMap;
  adaptiveHoleMap = new ADAPTIVE_HOLEMAP[maxtag];
  for(mi=0; mi<maxtag; mi++){
    adaptiveHoleMap[mi].meta.nlevel = 0;
    adaptiveHoleMap[mi].meta.elem_count = 0;
    adaptiveHoleMap[mi].meta.leaf_count = 0;
  }

  ADAPTIVE_HOLEMAP_COMPOSITE *adaptiveHoleMapCOMPOSITE;
  if(ncomposite){
    adaptiveHoleMapCOMPOSITE = new ADAPTIVE_HOLEMAP_COMPOSITE[maxtag];
    for(mi=0; mi<maxtag; mi++) adaptiveHoleMapCOMPOSITE[mi].meta.nlevel = 0;
  }

  /* ====================== */
  /* C: set existHole flags */
  /* ====================== */
  // determine holes for all bodies
  std::vector<uint8_t> existHole(maxtag);
  std::fill(existHole.begin(),existHole.end(),0);

  // fill array
  for(mbi=0; mbi<nblocks; mbi++){
    auto& mb = mblocks[mbi];
    existHole[mb->getMeshTag()-BASE] = (uint8_t) mb->getWallFlag();
  }
  MPI_Allreduce(MPI_IN_PLACE,existHole.data(),maxtag,MPI_UINT8_T,MPI_MAX,scomm);

  // set wall flags for all hole maps
  for(mi=0; mi<maxtag; mi++) adaptiveHoleMap[mi].existWall = existHole[mi];

  /* =============================== */
  /* D: construct Adaptive Hole Maps */
  /* =============================== */
  { // using encapsulation for clean up
    /* local variables */
    ADAPTIVE_HOLEMAP_OCTANT AHMO[nblocks];
    int i,j,l,c,n;

    for(mbi=0; mbi<nblocks; mbi++){
      ADAPTIVE_HOLEMAP_OCTANT &AHMOLocal = AHMO[mbi];
      auto& mb = mblocks[mbi];

      int adaptMap,adaptMapLocal;
      int nrefine;
      int existWallFlag;
      int meshtag;

      std::vector<char> refineFlag;
      std::vector<uint8_t> existWall;
      std::vector<uint8_t> existOuter;

      double bboxLocal[6];
      double bboxGlobal[6];
      double ds[3],dsmax,dsbox;

      /* =============================================== */
      /* Step 1: build global bounding box for this body */
      /* =============================================== */
      // get the local bounding box
      mb->getWallBounds(&meshtag,&existWallFlag,bboxLocal);

      // set hole for this body
      AHMOLocal.existWall = existHole[meshtag-BASE];

      // initialize global bounding box data
      for(i=0; i<3; i++) bboxGlobal[i]  = BIGVALUE;
      for(i=0; i<3; i++) bboxGlobal[3+i]=-BIGVALUE;

      // get the global bounding box info for this body (note the communicator)
      MPI_Allreduce(mb->bboxLocalAHM,bboxGlobal,3,MPI_DOUBLE,MPI_MIN,mb->blockcomm);
      MPI_Allreduce(&(mb->bboxLocalAHM[3]),&(bboxGlobal[3]),3,MPI_DOUBLE,MPI_MAX,mb->blockcomm);

      /* ======================================================= */
      /* Step 2: initialize and build adaptive map for this body */
      /* ======================================================= */
      // set extents and initialize map for this body
      if(AHMOLocal.existWall){
        for(j=0; j<3; j++){
          AHMOLocal.extents_lo[j] = bboxGlobal[j];
          AHMOLocal.extents_hi[j] = bboxGlobal[j+3];
          ds[j] = fabs(AHMOLocal.extents_hi[j] - AHMOLocal.extents_lo[j]);
        }

        // add buffer to extents
        dsmax = std::max(ds[0],ds[1]);
        dsmax = std::max(dsmax,ds[2]);
        dsbox = dsmax*0.01;
        for(j=0; j<3; j++) AHMOLocal.extents_lo[j] -= (dsbox);
        for(j=0; j<3; j++) AHMOLocal.extents_hi[j] += (dsbox);

        // initialize map for this body
        AHMOLocal.nlevel = 1;
        AHMOLocal.levels[0].elem_count = 0;
      }

      // initialize level 0, check outer boundary and set adapt flag
      adaptMapLocal = adaptMap = 0;
      if(AHMOLocal.existWall){
        level_octant_t *lvl = &AHMOLocal.levels[0];

        // allocate level 0
        lvl->level_id = 0;
        lvl->elem_count = 1;
        lvl->octants.resize(1);

        existWall.resize(1);
        existOuter.resize(1);

        // initialize
        existWall[0] = 1;
        existOuter[0] = 0;

        // fill in level 0 info
        lvl->octants[0].x = 0;
        lvl->octants[0].y = 0;
        lvl->octants[0].z = 0;
        lvl->octants[0].id = 1;
        lvl->octants[0].filltype = WALL_SB;
        lvl->octants[0].refined = 0;

        // set neighbor octants to NULL since Level 0
        for(n=0; n<6; n++) lvl->octants[0].nhbr[n] = NULL;

        // check if outer boundary exists
        mb->markBoundaryMapSurface(OUTERNODETYPE,
                                   AHMOLocal.extents_lo,
                                   AHMOLocal.extents_hi,
                                  &AHMOLocal.levels[0],
                                   NULL,
                                   existOuter.data());

        // check if initial octant contains both boundary types
        if(existOuter[0]){
          lvl->octants[0].refined = 1;
          adaptMapLocal = 1;
        }
      }
      // check this block's map for adaption (note the communicator)
      MPI_Allreduce(&adaptMapLocal,&adaptMap,1,MPI_INT,MPI_MAX,mb->blockcomm);

      // set all rank Outer flags
      if(adaptMap) AHMOLocal.levels[0].octants[0].refined = existOuter[0] = existWall[0] = 1;

      // recursively refine adaptive map until no intersection of wall and outer octants
      level_id = 0;
      while(adaptMap){
        level_octant_t *lvl = &AHMOLocal.levels[level_id];

        // update level count
        AHMOLocal.nlevel++;

        // reset adapt flag
        adaptMap = 0;

        // count number of octants to refine, set refined field before call octant_children
        refineFlag.resize(lvl->elem_count);

        nrefine = 0;
        for(i=0; i<lvl->elem_count; i++){
          refineFlag[i] = existWall[i] && existOuter[i];
          if(refineFlag[i]) nrefine++;
        }
        int nchildren = OCTANT_CHILDREN*nrefine;

        // allocate and fill next level
        level_id++;
        level_octant_t *new_lvl = &AHMOLocal.levels[level_id];

        // allocate new level
        new_lvl->level_id = level_id;
        new_lvl->elem_count = nchildren;
        new_lvl->octants.resize(nchildren);

        // zero arrays
        existWall.resize(nchildren);
        existOuter.resize(nchildren);
        std::fill(existWall.begin(),existWall.end(),0);
        std::fill(existOuter.begin(),existOuter.end(),0);

        // fill in new level children octant info
        nrefine = 0;
        for(i=0; i<lvl->elem_count; i++){
          if(refineFlag[i]){
            uint32_t newidx = OCTANT_CHILDREN*nrefine;

            // set Morton code and filltype for new octants
            octant_children(level_id,newidx,
                           &lvl->octants[i],
                           &new_lvl->octants[newidx+0],
                           &new_lvl->octants[newidx+1],
                           &new_lvl->octants[newidx+2],
                           &new_lvl->octants[newidx+3],
                           &new_lvl->octants[newidx+4],
                           &new_lvl->octants[newidx+5],
                           &new_lvl->octants[newidx+6],
                           &new_lvl->octants[newidx+7]);
            nrefine++;
          }
        }

        // fill in children neighbors AFTER ALL new level octants are formed
        nrefine = 0;
        for(i=0; i<lvl->elem_count; i++){
          if(refineFlag[i]){
            int newidx = OCTANT_CHILDREN*nrefine;

            // set neighbors for new octants
            octant_children_neighbors(&lvl->octants[i],
                                      &new_lvl->octants[newidx+0],
                                      &new_lvl->octants[newidx+1],
                                      &new_lvl->octants[newidx+2],
                                      &new_lvl->octants[newidx+3],
                                      &new_lvl->octants[newidx+4],
                                      &new_lvl->octants[newidx+5],
                                      &new_lvl->octants[newidx+6],
                                      &new_lvl->octants[newidx+7]);
            nrefine++;
          }
        }

        // check wall boundaries
        mb->markBoundaryMapSurface(WALLNODETYPE,
                                   AHMOLocal.extents_lo,
                                   AHMOLocal.extents_hi,
                                  &AHMOLocal.levels[level_id],
                                   NULL,
                                   existWall.data());

        // inform all mesh-block processes with this body tag of the octants flags (note the communicator)
        MPI_Allreduce(MPI_IN_PLACE,existWall.data(),nchildren,MPI_UINT8_T,MPI_MAX,mb->blockcomm);

        // check outer boundaries: check only octants with wall tagged previously
        mb->markBoundaryMapSurface(OUTERNODETYPE,
                                   AHMOLocal.extents_lo,
                                   AHMOLocal.extents_hi,
                                  &AHMOLocal.levels[level_id],
                                   existWall.data(),
                                   existOuter.data());

        // inform all mesh-block processes with this body tag of the octants flags (note the communicator)
        MPI_Allreduce(MPI_IN_PLACE,existOuter.data(),nchildren,MPI_UINT8_T,MPI_MAX,mb->blockcomm);

        // update filltype to wall if touching wall; check both boundary types
        for(i=0; i<new_lvl->elem_count; i++){
          if(existWall[i]){
            new_lvl->octants[i].filltype = WALL_SB;

            if(existOuter[i]){
              new_lvl->octants[i].refined = 1;
              adaptMap = 1;
            }
          }
        }
        if(AHMOLocal.nlevel >= OCTANT_MAXLEVEL && adaptMap == 1) {
          printf("[tioga] ERROR Constructing the Adaptive Hole Map! Max level reached!\n");
          exit(0);
        }
      }
      // adaption complete: all ranks of this mesh-block are informed of global octree
    }

    /* ======================================================== */
    /* Step 3: mark abutting/composite bodies holes in hole map */
    /* ======================================================== */
    for(int cb=0; cb<ncomposite; cb++){
      CompositeBody &Composite = compositeBody[cb];
      int nbodies = Composite.bodyids.size();

      for(i=0; i<nbodies; i++){
        int bodyi = Composite.bodyids[i]-BASE;
        ADAPTIVE_HOLEMAP_OCTANT &AHMOLocal = AHMO[bodyi];

        // check if this rank contains this body
        char rankContainsBody = 0;
        for(mbi=0; mbi<nblocks; mbi++){
          auto& mb = mblocks[mbi];
          int meshtag = mb->getMeshTag();
          if(meshtag-BASE == bodyi) {rankContainsBody = 1; break;}
        }

        meshblockCompInfo &MBC = meshblockComposite[bodyi];
        ADAPTIVE_HOLEMAP_COMPOSITE &AHMC = adaptiveHoleMapCOMPOSITE[bodyi];
        ahm_meta_minimal_t &meta = AHMC.meta;

        // initialize data for all ranks
        meta.nlevel = 0;

        /* --------------------------- */
        /* copy data for own mesh body */
        /* --------------------------- */
        if(MBC.comm != MPI_COMM_NULL){
          if(MBC.id == MBC.masterID){
            if(AHMOLocal.existWall){
              // 1. copy nlevel and extents
              meta.nlevel = AHMOLocal.nlevel;
              memcpy(meta.extents_lo,AHMOLocal.extents_lo,3*sizeof(double));
              memcpy(meta.extents_hi,AHMOLocal.extents_hi,3*sizeof(double));

              // 2. loop each level in adaptive hole map and assemble octants
              uint64_t elem_count = 0;
              for(level_id=0; level_id<meta.nlevel; level_id++){
                level_octant_t *lvl = &AHMOLocal.levels[level_id];
                level_octant_coordinate_t *elvl = &AHMC.levels[level_id];

                // set level info
                elvl->level_id = level_id;
                elvl->elem_count = lvl->elem_count;
                elvl->octants.resize(elvl->elem_count);

                // count octants
                elem_count += lvl->elem_count;

                // fill octant coordinate data
                for(j=0; j<elvl->elem_count; j++){
                  elvl->octants[j].x = lvl->octants[j].x;
                  elvl->octants[j].y = lvl->octants[j].y;
                  elvl->octants[j].z = lvl->octants[j].z;
                }
              }
              // 3. set element count
              meta.elem_count = elem_count;
            }
          }

          // 2. communicate minimal adaptive hole map to composite ranks
          // a. inform of meta data for map: packed into contiguous buffer
          MPI_Bcast(&meta,
                    sizeof(ahm_meta_minimal_t),
                    MPI_BYTE,
                    MBC.masterID,
                    MBC.comm);

          for(level_id=0; level_id<meta.nlevel; level_id++){
            level_octant_coordinate_t *level = &AHMC.levels[level_id];
            level->level_id = level_id;

            // communicate number of leaf octants on level
            MPI_Bcast(&(level->elem_count),
                      1,MPI_UINT32_T,
                      MBC.masterID,
                      MBC.comm);

            // allocate leaf octant data for level data (if not master rank)
            if(MBC.masterID != MBC.id) level->octants.resize(level->elem_count);
          }

          // d. loop each level in adaptive hole map: communicate octant coordinates
          //    - separate from previous loop to avoid MPI stalling
          for(level_id=0; level_id<meta.nlevel; level_id++){
            level_octant_coordinate_t *level = &AHMC.levels[level_id];

            // communicate leaf octant coordinate data on level
            MPI_Bcast(level->octants.data(),
                      level->elem_count*sizeof(octant_coordinates_t),
                      MPI_BYTE,
                      MBC.masterID,
                      MBC.comm);
          }

          // e. loop each level: check composite map wall boundaries for all mesh blocks
          for(level_id=0; level_id<meta.nlevel; level_id++){
            level_octant_coordinate_t *level = &AHMC.levels[level_id];
            std::vector<uint8_t> existWall(level->elem_count,0);

            // check composite map wall boundaries for all mesh blocks
            for(mbi=0; mbi<nblocks; mbi++){
              auto& mb = mblocks[mbi];

              // check if mesh block belongs to composite body
              if(compositeBodyMap[cb][bodyi][mb->getMeshTag()-BASE]){
                // check wall boundaries on abutting composite meshes
                mb->markBoundaryAdaptiveMapSurfaceIntersect(WALLNODETYPE,
                                                            meta.extents_lo,
                                                            meta.extents_hi,
                                                            level->level_id,
                                                            level->elem_count,
                                                            level->octants.data(),
                                                            NULL,
                                                            existWall.data());
              }
            }

            // f. inform all mesh-block processes with this body tag of the octants flags (note the communicator)
            MPI_Allreduce(MPI_IN_PLACE,existWall.data(),level->elem_count,MPI_UINT8_T,MPI_MAX,MBC.comm);

            // g. update filltype to wall if touching wall
            if(rankContainsBody){
              level_octant_t *new_lvl = &AHMOLocal.levels[level_id];
              for(int ii=0; ii<new_lvl->elem_count; ii++){
                if(existWall[ii]) new_lvl->octants[ii].filltype = WALL_SB;
              }
            }
          }
          // end marking wall BC in hole map for composite rank bodies
        }
      }
    }

    /* ======================================================== */
    /* Step 4: flood fill adaptive hole map for this mesh-block */
    /*   NOTE: all octants touching the octree boundary are     */
    /*         already tagged as OUTSIDE_SB (octant_children)   */
    /*         or as WALL_SB.                                   */
    /* ======================================================== */
    for(mbi=0; mbi<nblocks; mbi++){
      ADAPTIVE_HOLEMAP_OCTANT &AHMOLocal = AHMO[mbi];

      if(AHMOLocal.existWall){
        for(l=0; l<AHMOLocal.nlevel; l++) floodfill_level(&AHMOLocal.levels[l]);
      }
    }

    /* ================================================== */
    /* Step 5: build hole map and inform all complement   */
    /*         ranks of each mesh-block adaptive hole map */
    /* ================================================== */
    /* -------------------------------------------------- */
    /* assemble local adaptive hole map for own mesh body */
    /* -------------------------------------------------- */
    for(mbi=0; mbi<nblocks; mbi++){
      ADAPTIVE_HOLEMAP_OCTANT &AHMOLocal = AHMO[mbi];
      auto& mb = mblocks[mbi];
      int meshtag = mb->getMeshTag();

      ADAPTIVE_HOLEMAP &AHMLocal = adaptiveHoleMap[meshtag-BASE];
      ahm_meta_t &meta = AHMLocal.meta;

      // initialize meta
      meta.nlevel = 0;
      meta.leaf_count = 0;

      /* --------------------------- */
      /* copy data for own mesh body */
      /* --------------------------- */
      if(AHMLocal.existWall){
        // 1. copy nlevel and extents
        meta.nlevel = AHMOLocal.nlevel;
        memcpy(meta.extents_lo,AHMOLocal.extents_lo,3*sizeof(double));
        memcpy(meta.extents_hi,AHMOLocal.extents_hi,3*sizeof(double));

        // 2. loop each level in adaptive hole map and assemble octants
        for(level_id=0; level_id<meta.nlevel; level_id++){
          level_octant_t *lvl = &AHMOLocal.levels[level_id];
          level_t *elvl = &AHMLocal.levels[level_id];

          // set level info
          elvl->level_id = level_id;
          elvl->elem_count = lvl->elem_count;
          elvl->octants.resize(elvl->elem_count);

          // fill octant data
          for(j=0; j<elvl->elem_count; j++){
            elvl->octants[j].x = lvl->octants[j].x;
            elvl->octants[j].y = lvl->octants[j].y;
            elvl->octants[j].z = lvl->octants[j].z;
            elvl->octants[j].filltype = lvl->octants[j].filltype;
            elvl->octants[j].leafflag =!lvl->octants[j].refined;
            if(lvl->octants[j].refined){
              for(c=0; c<OCTANT_CHILDREN; c++) elvl->octants[j].children[c] =
                                                lvl->octants[j].children[c]->id;
            }

            // update leaf counter
            meta.leaf_count += elvl->octants[j].leafflag;
          }

          // free level full octant data set
          lvl->octants.clear();
        }
      }
    }
  }

  /* ========================================================= */
  /* E. communicate adaptive hole map info to complement ranks */
  /* ========================================================= */
  for(int i=0; i<maxtag; i++){
    ADAPTIVE_HOLEMAP &AHME = adaptiveHoleMap[i];
    ahm_meta_t &meta = AHME.meta;

    meta.elem_count = 0;
    if(AHME.existWall){
      meshblockCompInfo &MBC = meshblockComplement[i];

      // complement + master ranks only involved
      if(MBC.comm != MPI_COMM_NULL){

        // 1. inform of meta data for map: packed into contiguous buffer
        MPI_Bcast(&AHME.meta,
                  sizeof(ahm_meta_t),
                  MPI_BYTE,
                  MBC.masterID,
                  MBC.comm);

        // 2. loop each level in adaptive hole map and assemble octants
        for(level_id=0; level_id<meta.nlevel; level_id++){
          level_t *elvl = &AHME.levels[level_id];

          // set level id
          elvl->level_id = level_id;

          // communicate number of leaf octants on level
          MPI_Bcast(&(elvl->elem_count),
                    1,MPI_UINT32_T,
                    MBC.masterID,
                    MBC.comm);

          // allocate leaf octant data for level data (if not master rank)
          if(MBC.masterID != MBC.id) elvl->octants.resize(elvl->elem_count);

          // communicate leaf octant data on level
          MPI_Bcast(elvl->octants.data(),
                    elvl->elem_count*sizeof(octant_t),
                    MPI_BYTE,
                    MBC.masterID,
                    MBC.comm);

          // count number of octants
          meta.elem_count += elvl->elem_count;
        }
      }
    } else {
      meta.leaf_count = 0;
    }
  }
  MPI_Barrier(scomm);

  // clean up memory
  if(ncomposite) delete [] adaptiveHoleMapCOMPOSITE;

  // set the global number of meshes to maxtag
  nmesh=maxtag;

  // output statistics
  if(myid==0){
    for(int i=0;i<nmesh;i++){
      if(adaptiveHoleMap[i].meta.nlevel > 0){
        fprintf(stdout,"  "
              "[tioga::performConnectivity::getAdaptiveHoleMap] "
              "Mesh Body %d: Levels %d, "
              "total octants %lu, total leafs %lu\n",i,
              adaptiveHoleMap[i].meta.nlevel,
              adaptiveHoleMap[i].meta.elem_count,
              adaptiveHoleMap[i].meta.leaf_count);
        fflush(stdout);
      }
    }
  }

  // output the hole maps
//  for(mbi=0; mbi<nblocks; mbi++){
//    auto& mb = mblocks[mbi];
//    mb->writeBCnodes(WALLNODETYPE,mb->getMeshTag()-BASE);
//  }
//  this->outputAdaptiveHoleMap();
}

/**
 * Output the hole map to a tecplot compatible file
*/
void tioga::outputHoleMap(void)
{
  int i,k;
  int nnodes,ncells;
  int ns1,ns2;
  int ii,jj,kk,m;
  FILE *fp;
  double ds[3];
  char intstring[12];
  char fname[80];

  for(i=0;i<nmesh;i++)
    if (holeMap[i].existWall)
       {
	 sprintf(intstring,"%d",100000+i+100*myid);
	 sprintf(fname,"holeMap%s.dat",&(intstring[1]));
	 fp=fopen(fname,"w");
	 fprintf(fp,"TITLE =\"Tioga output\"\n");
	 fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\"\n");
	 nnodes=(holeMap[i].nx[0]+1)*(holeMap[i].nx[1]+1)*(holeMap[i].nx[2]+1);
	 ncells=(holeMap[i].nx[0])*(holeMap[i].nx[1])*(holeMap[i].nx[2]);
	 fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEBLOCK\n",nnodes,ncells);
	 fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED)\n");
	 for(k=0;k<3;k++) ds[k]=(holeMap[i].extents[k+3]-holeMap[i].extents[k])/(holeMap[i].nx[k]);
	 //
	 for(kk=0;kk<holeMap[i].nx[2]+1;kk++)
	   for(jj=0;jj<holeMap[i].nx[1]+1;jj++)
	     for(ii=0;ii<holeMap[i].nx[0]+1;ii++)
	       fprintf(fp,"%.14e\n",ii*ds[0]);
	 for(kk=0;kk<holeMap[i].nx[2]+1;kk++)
	   for(jj=0;jj<holeMap[i].nx[1]+1;jj++)
	     for(ii=0;ii<holeMap[i].nx[0]+1;ii++)
	       fprintf(fp,"%.14e\n",jj*ds[1]);
	 for(kk=0;kk<holeMap[i].nx[2]+1;kk++)
	   for(jj=0;jj<holeMap[i].nx[1]+1;jj++)
	     for(ii=0;ii<holeMap[i].nx[0]+1;ii++)
	       fprintf(fp,"%.14e\n",kk*ds[2]);
	 m=0;
	 for(kk=0;kk<holeMap[i].nx[2];kk++)
	   for(jj=0;jj<holeMap[i].nx[1];jj++)
	     for(ii=0;ii<holeMap[i].nx[0];ii++)
	       {
		 fprintf(fp,"%f\n",(double)holeMap[i].sam[m]);
		 m++;
	       }

	 m=0;
         ns1=holeMap[i].nx[0]+1;
	 ns2=(holeMap[i].nx[1]+1)*ns1;
	 for(kk=0;kk<holeMap[i].nx[2];kk++)
	   for(jj=0;jj<holeMap[i].nx[1];jj++)
	     for(ii=0;ii<holeMap[i].nx[0];ii++)
	       {
		 m=kk*ns2+jj*ns1+ii+1;
		 fprintf(fp,"%d %d %d %d %d %d %d %d\n",m,m+1,m+1+ns1,m+ns1,
			 m+ns2,m+1+ns2,m+ns2+ns1+1,m+ns1+ns2);
	       }
       }
 fclose(fp);
}

void writePointsHeaderVolume(const char *filename){
  FILE *fp = fopen(filename,"w");
  fprintf(fp,"TITLE =\"Octree Volume Points\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"IBLANK\",\"LEVEL\",\"ID\"\n");
  fflush(fp);
  fclose(fp);
}

void writePointsVolume(FILE *fp,int level,int id,double *x,int npts1d,int type){
  int i,j,k;
  int pt;

  /* write zone header */
  fprintf(fp,"ZONE T=\"Volume_L%d_Q%d\", I=%d, J=%d, K=%d, F=POINT\n",
          level,id,npts1d,npts1d,npts1d);

  for (pt=0,k=0; k<npts1d; ++k) {
    for (j=0; j<npts1d; ++j) {
      for (i=0; i<npts1d; ++i) {
        double *co = &x[3*pt];
        fprintf(fp,"%f %f %f %d %d %d\n",co[0],co[1],co[2],type,level,id);
        pt++;
      }
    }
  }
  fflush(fp);
}

void writePointsVolumeBrick(FILE *fp,int level,int id,double *x,int npts1d,int type){
  int i,j,k;
  int pt;

  /* write zone header */
  fprintf(fp,"ZONE T=\"Volume_L%d_Q%d\",N=%d, E=%d ET=BRICK, F=FEBLOCK\n",level,id,8,1);
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL, 4=CELLCENTERED)\n");

  // x points
  for (pt=0,k=0; k < npts1d; ++k) {
    for (j=0; j < npts1d; ++j) {
      for (i = 0; i < npts1d; ++i) {
        double *co = &x[3*pt];
        fprintf(fp,"%f ",co[0]);
        pt++;
      }
    }
  }
  fprintf(fp,"\n");

  // y points
  for (pt=0,k=0; k < npts1d; ++k) {
    for (j=0; j < npts1d; ++j) {
      for (i = 0; i < npts1d; ++i) {
        double *co = &x[3*pt];
        fprintf(fp,"%f ",co[1]);
        pt++;
      }
    }
  }
  fprintf(fp,"\n");

  // z points
  for (pt=0,k=0; k < npts1d; ++k) {
    for (j=0; j < npts1d; ++j) {
      for (i = 0; i < npts1d; ++i) {
        double *co = &x[3*pt];
        fprintf(fp,"%f ",co[2]);
        pt++;
      }
    }
  }
  fprintf(fp,"\n");

  // data
  fprintf(fp,"%d\n",type);

  // connectivity: list rectangles nodes
  fprintf(fp,"1 2 4 3\n");
  fprintf(fp,"5 6 8 7\n");
  fprintf(fp,"1 5 7 3\n");
  fprintf(fp,"2 6 8 4\n");
  fprintf(fp,"1 2 6 5\n");
  fprintf(fp,"3 7 8 4\n");
  fflush(fp);
}

/**
 * Output the adaptive hole map to a tecplot compatible file
 */
void tioga::outputAdaptiveHoleMap(void){
  char filename[128];
  FILE *file;

  double pt_coords[3*8];
  double ds[3],dx[3];
  double x[2],y[2],z[2];
  int level,e,m;
  int i,j,k;
  int pt;

  static int ahm_step = 0;
  if(myid==0){
    for(m=0;m<nmesh;m++){
      if(adaptiveHoleMap[m].existWall){
        ahm_meta_t &meta = adaptiveHoleMap[m].meta;

        ds[0] = meta.extents_hi[0] - meta.extents_lo[0];
        ds[1] = meta.extents_hi[1] - meta.extents_lo[1];
        ds[2] = meta.extents_hi[2] - meta.extents_lo[2];

        sprintf(filename,"AHM.body%d.%d.tec",m,ahm_step++);
        writePointsHeaderVolume(filename);

        file = fopen(filename, "a");
        for(level = (meta.nlevel>1) ? 1:0; level < meta.nlevel; level++){
          level_t &L = adaptiveHoleMap[m].levels[level];
          const qcoord_t levelh = OCTANT_LEN(L.level_id);

          // octant length
          dx[0] = ds[0]*INT2DBL*levelh;
          dx[1] = ds[1]*INT2DBL*levelh;
          dx[2] = ds[2]*INT2DBL*levelh;

          for(e = 0; e < L.elem_count; e++){
            if(L.octants[e].leafflag==0) continue; // skip to next element

            // compute physical coordinates
            x[0] = meta.extents_lo[0] + ds[0]*INT2DBL*L.octants[e].x;
            y[0] = meta.extents_lo[1] + ds[1]*INT2DBL*L.octants[e].y;
            z[0] = meta.extents_lo[2] + ds[2]*INT2DBL*L.octants[e].z;
            x[1] = x[0] + dx[0];
            y[1] = y[0] + dx[1];
            z[1] = z[0] + dx[2];

            pt = 0;
            for(k=0;k<2;k++){
              for(j=0;j<2;j++){
                for(i=0;i<2;i++){
                  pt_coords[3*pt+0] = x[i];
                  pt_coords[3*pt+1] = y[j];
                  pt_coords[3*pt+2] = z[k];
                  pt++;
                }
              }
            }
            writePointsVolume(file,level,e,pt_coords,2,L.octants[e].filltype);
          }
        }
        fflush(file);
        fclose(file);
      }
    }
  }
}