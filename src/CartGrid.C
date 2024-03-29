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
#include "tioga_gpu.h"
#include "TiogaMeshInfo.h"
#include "codetypes.h"
#include "CartGrid.h"

#include <cassert>
#include <numeric>

CartGrid::~CartGrid()
{
  if (own_data_ptrs) {
    if (global_id) TIOGA_FREE(global_id);
    if (level_num) TIOGA_FREE(level_num);
    if (proc_id) TIOGA_FREE(proc_id);
    if (local_id) TIOGA_FREE(local_id);
    if (ilo) TIOGA_FREE(ilo);
    if (ihi) TIOGA_FREE(ihi);
    if (dims) TIOGA_FREE(dims);
    if (xlo) TIOGA_FREE(xlo);
    if (dx) TIOGA_FREE(dx);
  }

  if (own_amr_mesh_info && (m_info != nullptr)) {
    TIOGA_FREE_DEVICE(m_info->level.dptr);
    TIOGA_FREE_DEVICE(m_info->mpi_rank.dptr);
    TIOGA_FREE_DEVICE(m_info->local_id.dptr);
    TIOGA_FREE_DEVICE(m_info->ilow.dptr);
    TIOGA_FREE_DEVICE(m_info->ihigh.dptr);
    TIOGA_FREE_DEVICE(m_info->dims.dptr);
    TIOGA_FREE_DEVICE(m_info->xlo.dptr);
    TIOGA_FREE_DEVICE(m_info->dx.dptr);

    delete m_info;
  }

  if (lcount) TIOGA_FREE(lcount);
  if (dxlvl) TIOGA_FREE(dxlvl);

  if (m_info_device != nullptr) TIOGA_FREE_DEVICE(m_info_device);
};

void CartGrid::registerData(TIOGA::AMRMeshInfo* minfo)
{
  own_data_ptrs = false;
  own_amr_mesh_info = false;
  m_info = minfo;
  ngrids = minfo->ngrids_global;
  global_id = nullptr; // unused
  level_num = minfo->level.hptr;
  proc_id = minfo->mpi_rank.hptr;
  local_id = minfo->local_id.hptr;
  ilo = minfo->ilow.hptr;
  ihi = minfo->ihigh.hptr;
  dims = minfo->dims.hptr;
  xlo = minfo->xlo.hptr;
  dx = minfo->dx.hptr;
  nf = minfo->num_ghost;

  if (m_info_device == nullptr) {
    m_info_device = TIOGA::gpu::allocate_on_device<TIOGA::AMRMeshInfo>(
      sizeof(TIOGA::AMRMeshInfo));
  }
  TIOGA::gpu::copy_to_device(m_info_device, m_info, sizeof(TIOGA::AMRMeshInfo));
}

void CartGrid::registerData(int nfin,int *idata,double *rdata,int ngridsin)
{
  int i,i3,i6,iloc,n;
  FILE *fp;
  ngrids = ngridsin;
  global_id=(int *) malloc(sizeof(int)*ngrids);
  level_num=(int *) malloc(sizeof(int)*ngrids);
  proc_id=(int *) malloc(sizeof(int)*ngrids);
  ilo=(int *) malloc(sizeof(int)*3*ngrids);
  ihi=(int *) malloc(sizeof(int)*3*ngrids);
  xlo=(double *) malloc(sizeof(double)*3*ngrids);
  dx=(double *) malloc(sizeof(double)*3*ngrids);
  local_id=(int *)malloc(sizeof(int)*ngrids);
  dims=(int *)malloc(sizeof(dims)*3*ngrids);
  nf=nfin;
  if (myid==0) fp=fopen("cartGrid.dat","w");
  for(i=0;i<ngrids;i++)
    {
      i3=3*i;
      i6=2*i3;
      iloc=10*i;

      global_id[i]=idata[iloc];
      level_num[i]=idata[iloc+1];
      proc_id[i]=idata[iloc+2];
      local_id[i]=idata[iloc+3];
      for(n=0;n<3;n++)
	{
	  ilo[i3+n]=idata[iloc+4+n];
	  ihi[i3+n]=idata[iloc+7+n];
	  dims[i3+n]=ihi[i3+n]-ilo[i3+n]+1;
	}
      xlo[i3]=rdata[i6];
      xlo[i3+1]=rdata[i6+1];
      xlo[i3+2]=rdata[i6+2];
      dx[i3]=rdata[i6+3];
      dx[i3+1]=rdata[i6+4];
      dx[i3+2]=rdata[i6+5];
      if (myid==0) 
        fprintf(fp,"%d %d %d %d %f %f %f\n",global_id[i],level_num[i],proc_id[i],
                                   local_id[i],dx[i3],dx[i3+1],dx[i3+2]);
    }
   if (myid==0) fclose(fp);

   // Create AMRMeshInfo object so that it can be accessed on device in future
   // create_mesh_info();
};

//
// Bare bone preprocessor now
// willl add more support data structures
// to promote efficient search once the concept
// works
//
void CartGrid::preprocess(void)
{
  int i,n;
  //
  // find the global minimum coord location
  //
  xlosup[0]=xlosup[1]=xlosup[2]=BIGVALUE;
  maxlevel=-1;
  for (i=0;i<ngrids;i++)
    {
      for(n=0;n<3;n++)
	xlosup[n]=((xlosup[n] <= xlo[3*i+n]) ? xlosup[n]: xlo[3*i+n]);
      maxlevel=((maxlevel >= level_num[i]) ? maxlevel: level_num[i]);
    }
    maxlevel++;
  lcount=(int *)malloc(sizeof(int)*maxlevel);
  dxlvl=(double *)malloc(sizeof(double)*3*maxlevel);
  for(i=0;i<maxlevel;i++) lcount[i]=0;
  for(i=0;i<ngrids;i++)
    {
      lcount[level_num[i]]++;
      for(n=0;n<3;n++)
	dxlvl[3*level_num[i]+n]=dx[3*i+n];
    }
}
//
// Basic search routine now
// will improve efficiency once it works
//
void CartGrid::search(double *x,int *donorid,int npts)
{
  int i,j,k,l,n,il[3];
  bool flag;
  int dcount;
  dcount=0;
  for(i=0;i<npts;i++)
    {
      flag=0;
      donorid[i]=-1;
      for(l=maxlevel-1;l>=0 && flag==0;l--)
	{
	  for(n=0;n<3;n++)
	    il[n]=floor((x[3*i+n]-xlosup[n])/dxlvl[3*l+n]);
	  for(j=0;j<ngrids && flag==0;j++)
	    {
	      if (level_num[j]==l) 
		{
		  flag=1;
		  // for(n=0;n<3;n++) flag=flag && (x[3*i+n] >=xlo[3*j+n]);
		  // for(n=0;n<3;n++) flag=flag && (x[3*i+n] <=xlo[3*j+n]+
		  //   			 dx[3*j+n]*(dims[3*j+n]));
		  //for(n=0;n<3;n++) flag = flag && (il[n] >=ilo[3*j+n]);
		  //for(n=0;n<3;n++) flag = flag && (il[n] <=ihi[3*j+n]);
          for(n=0;n<3;n++) flag=flag && ((x[3*i+n]-xlo[3*j+n]) > -TOL);
          for(n=0;n<3;n++) flag=flag && ((x[3*i+n]- (xlo[3*j+n]+
                                                     dx[3*j+n]*(dims[3*j+n]))) < TOL);
		  if (flag) { 
		    dcount++; 
		    donorid[i]=j; 
		}
	    }
	}
    }
    if (myid==2 && abs(x[3*i]-0.739573) < 1e-5 && abs(x[3*i+1]+0.259310) < 1e-5 &&
        abs(x[3*i+2]+0.639614) < 1e-5) {
     printf("%d %d %f %f %f %d\n",myid,i,x[3*i],x[3*i+1],x[3*i+2],donorid[i]);
   }
    if (donorid[i]==-1) printf("%d %f %f %f\n",myid,x[3*i],x[3*i+1],x[3*i+2]);
  }
 //printf("CartGrid::search Processor %d located %d of %d points\n",myid,dcount,npts);
}

namespace {

template<typename T>
inline void create_view(TIOGA::TiogaView<T> tv, T* sptr, int sz)
{
  tv.sz = sz;
  tv.hptr = sptr;
  tv.dptr = TIOGA::gpu::push_to_device<T>(tv.hptr, sizeof(T) * sz);
}

}

/** Create AMRMeshInfo objects that can be accessed on host and device
 *
 *  This method serves two purposes: 1. It ensures that a valid AMRMeshInfo
 *  object exists even when using the legacy TIOGA API; and, 2. For simulations
 *  where AMR solver does not span all MPI processes visible to TIOGA, this
 *  creates a valid AMRMeshInfo object so that the unstructured block can query
 *  patch information to perform searches. See tioga::preprocess_amr_data for
 *  more information.
 */
void CartGrid::create_mesh_info()
{
  assert(proc_id != nullptr);

  if (m_info == nullptr) m_info = new TIOGA::AMRMeshInfo;
  m_info->ngrids_global = ngrids;
  m_info->num_ghost = nf;
  create_view(m_info->level, level_num, ngrids);
  create_view(m_info->mpi_rank, proc_id, ngrids);
  create_view(m_info->local_id, local_id, ngrids);
  create_view(m_info->ilow, ilo, ngrids * 3);
  create_view(m_info->ihigh, ihi, ngrids * 3);
  create_view(m_info->dims, dims, ngrids * 3);
  create_view(m_info->xlo, xlo, ngrids * 3);
  create_view(m_info->dx, dx, ngrids * 3);

  int iproc = myid;
  int nplocal = std::accumulate(
    proc_id, proc_id + ngrids, 0,
    [iproc](int x, int y) -> int { return x + ((iproc == y) ? 1 : 0); });

  m_info->ngrids_local = nplocal;

  if (m_info_device == nullptr) {
    m_info_device = TIOGA::gpu::allocate_on_device<TIOGA::AMRMeshInfo>(
      sizeof(TIOGA::AMRMeshInfo));
  }
  TIOGA::gpu::copy_to_device(m_info_device, m_info, sizeof(TIOGA::AMRMeshInfo));

  own_data_ptrs = true;
  own_amr_mesh_info = true;
}
