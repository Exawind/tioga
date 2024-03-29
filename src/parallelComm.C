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
#include "codetypes.h"
#include "mpi.h"
#include <algorithm>
#include <numeric>
#include "parallelComm.h"
#define REAL double

void parallelComm::sendRecvPacketsAll(PACKET *sndPack, PACKET *rcvPack)
{
  int i;
  int *sint,*sreal,*rint,*rreal;
  //
  sint=(int *)malloc(sizeof(int)*numprocs);
  sreal=(int *) malloc(sizeof(int)*numprocs);
  rint=(int *)malloc(sizeof(int)*numprocs);
  rreal=(int *) malloc(sizeof(int)*numprocs);
  //
  for(i=0;i<numprocs;i++){
    sint[i]=sndPack[i].nints;			
    sreal[i]=sndPack[i].nreals;
  }
  //
  MPI_Alltoall(sint,1,MPI_INT,rint,1,MPI_INT,scomm);
  MPI_Alltoall(sreal,1,MPI_INT,rreal,1,MPI_INT,scomm);
  //
  for(i=0;i<numprocs;i++) {
    rcvPack[i].nints=rint[i];
    rcvPack[i].nreals=rreal[i];
  }
  //

  int all_snd_nints = std::accumulate(sint, sint + numprocs, 0);
  int all_rcv_nints = std::accumulate(rint, rint + numprocs, 0);
  int *all_snd_intData, *all_rcv_intData;
  all_snd_intData=(int *) malloc(sizeof(int)*all_snd_nints);
  all_rcv_intData=(int *) malloc(sizeof(int)*all_rcv_nints);
  std::vector<int> snd_int_displs(numprocs+1, 0);
  std::vector<int> rcv_int_displs(numprocs+1, 0);
  for (int i=1; i <= numprocs; i++) {
    snd_int_displs[i] = snd_int_displs[i-1] + sint[i-1];
    rcv_int_displs[i] = rcv_int_displs[i-1] + rint[i-1];
  }
  for (int i=0; i < numprocs; i++) {
    int displ = snd_int_displs[i];
    for(int j=0; j < sint[i]; j++){
      all_snd_intData[displ+j] = sndPack[i].intData[j];
    }
  }
  MPI_Request int_request;
  MPI_Ialltoallv(all_snd_intData,
                sint,
                snd_int_displs.data(),
                MPI_INT,
                all_rcv_intData,
                rint,
                rcv_int_displs.data(),
                MPI_INT,
                scomm,
                &int_request);

  int all_snd_nreals = std::accumulate(sreal, sreal + numprocs, 0);
  int all_rcv_nreals = std::accumulate(rreal, rreal + numprocs, 0);
  REAL *all_snd_realData, *all_rcv_realData;
  all_snd_realData=(REAL *) malloc(sizeof(REAL)*all_snd_nreals);
  all_rcv_realData=(REAL *) malloc(sizeof(REAL)*all_rcv_nreals);
  std::vector<int> snd_real_displs(numprocs+1, 0);
  std::vector<int> rcv_real_displs(numprocs+1, 0);
  for (int i=1; i <= numprocs; i++) {
    snd_real_displs[i] = snd_real_displs[i-1] + sreal[i-1];
    rcv_real_displs[i] = rcv_real_displs[i-1] + rreal[i-1];
  }
  for (int i=0; i < numprocs; i++) {
    int displ = snd_real_displs[i];
    for(int j=0; j < sreal[i]; j++){
      all_snd_realData[displ+j] = sndPack[i].realData[j];
    }
  }
  MPI_Request real_request;
  MPI_Ialltoallv(all_snd_realData,
                sreal,
                snd_real_displs.data(),
                MPI_DOUBLE,
                all_rcv_realData,
                rreal,
                rcv_real_displs.data(),
                MPI_DOUBLE,
                scomm,
                &real_request);

  MPI_Wait(&int_request, MPI_STATUS_IGNORE);
  for(i=0;i<numprocs;i++){
    if (rcvPack[i].nints > 0) {
      rcvPack[i].intData=(int *) malloc(sizeof(int)*rcvPack[i].nints);
    }
    if (rcvPack[i].nreals > 0) {
      rcvPack[i].realData=(REAL *) malloc(sizeof(REAL)*rcvPack[i].nreals);
    }
  }

  MPI_Wait(&real_request, MPI_STATUS_IGNORE);
  for (int i=0; i < numprocs; i++) {
    int displ = rcv_int_displs[i];
    for(int j=0; j < rint[i]; j++){
      rcvPack[i].intData[j] = all_rcv_intData[displ+j];
    }
  }
  for (int i=0; i < numprocs; i++) {
    int displ = rcv_real_displs[i];
    for(int j=0; j < rreal[i]; j++){
      rcvPack[i].realData[j] = all_rcv_realData[displ+j];
    }
  }

  TIOGA_FREE(all_snd_intData);
  TIOGA_FREE(all_rcv_intData);
  TIOGA_FREE(all_snd_realData);
  TIOGA_FREE(all_rcv_realData);
  TIOGA_FREE(sint);
  TIOGA_FREE(sreal);
  TIOGA_FREE(rint);
  TIOGA_FREE(rreal);
}

// Old version that does not use ialltoallv
// void parallelComm::sendRecvPacketsAll(PACKET *sndPack, PACKET *rcvPack)
// {
//   int i;
//   int *sint,*sreal,*rint,*rreal;
//   int tag,irnum;
//   MPI_Request *request;
//   MPI_Status *status;
//   //
//   sint=(int *)malloc(sizeof(int)*numprocs);
//   sreal=(int *) malloc(sizeof(int)*numprocs);
//   rint=(int *)malloc(sizeof(int)*numprocs);
//   rreal=(int *) malloc(sizeof(int)*numprocs);
//   request=(MPI_Request *) malloc(sizeof(MPI_Request)*4*numprocs);
//   status=(MPI_Status *) malloc(sizeof(MPI_Status)*4*numprocs);
//   //
//   for(i=0;i<numprocs;i++){
//     sint[i]=sndPack[i].nints;			
//     sreal[i]=sndPack[i].nreals;
//   }
//   //
//   MPI_Alltoall(sint,1,MPI_INT,rint,1,MPI_INT,scomm);
//   MPI_Alltoall(sreal,1,MPI_INT,rreal,1,MPI_INT,scomm);
//   //
//   for(i=0;i<numprocs;i++) {
//     rcvPack[i].nints=rint[i];
//     rcvPack[i].nreals=rreal[i];
//   }
//   //
//   irnum=0;
//   for(i=0;i<numprocs;i++)
//     {
//       if (rcvPack[i].nints > 0) {
// 	tag=1;
// 	rcvPack[i].intData=(int *) malloc(sizeof(int)*rcvPack[i].nints);
// 	MPI_Irecv(rcvPack[i].intData,rcvPack[i].nints,
// 		  MPI_INT,i,
// 		  tag,scomm,&request[irnum++]);
//       }
//       if (rcvPack[i].nreals > 0) {
// 	tag=2;
// 	rcvPack[i].realData=(REAL *) malloc(sizeof(REAL)*rcvPack[i].nreals);
// 	MPI_Irecv(rcvPack[i].realData,rcvPack[i].nreals,
// 		  MPI_DOUBLE,i,
// 		  tag,scomm,&request[irnum++]);
//       }
//     }
//   for(i=0;i<numprocs;i++)
//     {
//       if (sndPack[i].nints > 0){
// 	tag=1;
// 	MPI_Isend(sndPack[i].intData,sndPack[i].nints,
// 		  MPI_INT,i,
// 		  tag,scomm,&request[irnum++]);
//       }
//       if (sndPack[i].nreals > 0){
// 	tag=2;
// 	MPI_Isend(sndPack[i].realData,sndPack[i].nreals,
// 		  MPI_DOUBLE,i,
// 		  tag,scomm,&request[irnum++]);
//       }
//     }
//   MPI_Waitall(irnum,request,status);
  
//   TIOGA_FREE(sint);
//   TIOGA_FREE(sreal);
//   TIOGA_FREE(rint);
//   TIOGA_FREE(rreal);
//   TIOGA_FREE(request);
//   TIOGA_FREE(status);
// }

void parallelComm::sendRecvPackets(PACKET *sndPack,PACKET *rcvPack)
{
  int i;
  int *scount,*rcount;
  int tag,irnum;
  MPI_Request *request;
  MPI_Status *status;
  //
  scount=(int *)malloc(2*sizeof(int)*nsend);
  rcount=(int *) malloc(2*sizeof(int)*nrecv);
  request=(MPI_Request *) malloc(sizeof(MPI_Request)*2*(nsend+nrecv));
  status=(MPI_Status *) malloc(sizeof(MPI_Status)*2*(nsend+nrecv));
  //
  for(i=0;i<nsend;i++){
    scount[2*i]=sndPack[i].nints;			
    scount[2*i+1]=sndPack[i].nreals;
  }
  //
  irnum=0;
  tag=1;
  //
  for(i=0;i<nrecv;i++)
    MPI_Irecv(&(rcount[2*i]),2,MPI_INT,rcvMap[i],tag,scomm,&request[irnum++]);
  //
  for(i=0;i<nsend;i++)
    MPI_Isend(&(scount[2*i]),2,MPI_INT,sndMap[i],tag,scomm,&request[irnum++]);
  //
  MPI_Waitall(irnum,request,status);
  for(i=0;i<nrecv;i++)
    {
      rcvPack[i].nints=rcount[2*i];
      rcvPack[i].nreals=rcount[2*i+1];
    }
  //
  irnum=0;
  for(i=0;i<nrecv;i++)
    {
      if (rcvPack[i].nints > 0) {
	tag=1;
	rcvPack[i].intData=(int *) malloc(sizeof(int)*rcvPack[i].nints);
	MPI_Irecv(rcvPack[i].intData,rcvPack[i].nints,
		  MPI_INT,rcvMap[i],
		  tag,scomm,&request[irnum++]);
      }
      if (rcvPack[i].nreals > 0) {
	tag=2;
	rcvPack[i].realData=(REAL *) malloc(sizeof(REAL)*rcvPack[i].nreals);
	MPI_Irecv(rcvPack[i].realData,rcvPack[i].nreals,
		  MPI_DOUBLE,rcvMap[i],
		  tag,scomm,&request[irnum++]);
      }
    }
  //
  for(i=0;i<nsend;i++)
    {
      if (sndPack[i].nints > 0){
	tag=1;
	MPI_Isend(sndPack[i].intData,sndPack[i].nints,
		  MPI_INT,sndMap[i],
		  tag,scomm,&request[irnum++]);
      }
      if (sndPack[i].nreals > 0){
	tag=2;
	MPI_Isend(sndPack[i].realData,sndPack[i].nreals,
		  MPI_DOUBLE,sndMap[i],
		  tag,scomm,&request[irnum++]);
      }
    }
  MPI_Waitall(irnum,request,status);
  //
  TIOGA_FREE(scount);
  TIOGA_FREE(rcount);
  TIOGA_FREE(request);
  TIOGA_FREE(status);
}

void parallelComm::sendRecvPacketsCheck(PACKET *sndPack,PACKET *rcvPack)
{
  int i;
  int *scount,*rcount;
  int tag,irnum;
  MPI_Request *request;
  MPI_Status *status;
  //
  scount=(int *)malloc(2*sizeof(int)*nsend);
  rcount=(int *) malloc(2*sizeof(int)*nrecv);
  request=(MPI_Request *) malloc(sizeof(MPI_Request)*2*(nsend+nrecv));
  status=(MPI_Status *) malloc(sizeof(MPI_Status)*2*(nsend+nrecv));
  //
  for(i=0;i<nsend;i++){
    scount[2*i]=sndPack[i].nints;			
    scount[2*i+1]=sndPack[i].nreals;
  }
  //
  irnum=0;
  tag=1;
  //
  for(i=0;i<nrecv;i++)
    MPI_Irecv(&(rcount[2*i]),2,MPI_INT,rcvMap[i],tag,scomm,&request[irnum++]);
  //
  for(i=0;i<nsend;i++)
    MPI_Isend(&(scount[2*i]),2,MPI_INT,sndMap[i],tag,scomm,&request[irnum++]);
  //
  MPI_Waitall(irnum,request,status);

  for(i=0;i<nrecv;i++)
    {
      rcvPack[i].nints=rcount[2*i];
      rcvPack[i].nreals=rcount[2*i+1];
    }

  //for(i=0;i<nsend;i++)
  //  {
  //    printf("%d sending %d to %d\n",myid,sndPack[i].nints,sndMap[i]);
  //  }	
  //for(i=0;i<nrecv;i++)
  //  {
  //   printf("%d receiving %d from %d\n",myid,rcvPack[i].nints,rcvMap[i]);
  //  }
  //
  irnum=0;
  for(i=0;i<nrecv;i++)
    {
      if (rcvPack[i].nints > 0) {
	tag=1;
	rcvPack[i].intData=(int *) malloc(sizeof(int)*rcvPack[i].nints);
	MPI_Irecv(rcvPack[i].intData,rcvPack[i].nints,
		  MPI_INT,rcvMap[i],
		  tag,scomm,&request[irnum++]);
      }
      if (rcvPack[i].nreals > 0 ) {
	tag=2;
	rcvPack[i].realData=(REAL *) malloc(sizeof(REAL)*rcvPack[i].nreals);
	MPI_Irecv(rcvPack[i].realData,rcvPack[i].nreals,
		  MPI_DOUBLE,rcvMap[i],
		  tag,scomm,&request[irnum++]);
      }
    }
  //
  for(i=0;i<nsend;i++)
    {
      if (sndPack[i].nints > 0){
	tag=1;
	MPI_Isend(sndPack[i].intData,sndPack[i].nints,
		  MPI_INT,sndMap[i],
		  tag,scomm,&request[irnum++]);
      }
      if (sndPack[i].nreals > 0){
	tag=2;
	MPI_Isend(sndPack[i].realData,sndPack[i].nreals,
		  MPI_DOUBLE,sndMap[i],
		  tag,scomm,&request[irnum++]);
      }
    }
  MPI_Waitall(irnum,request,status);
  //
  TIOGA_FREE(scount);
  TIOGA_FREE(rcount);
  TIOGA_FREE(request);
  TIOGA_FREE(status);
}

void parallelComm::setMap(int ns,int nr, int *snd,int *rcv)
{
  int i;
  //
  if (sndMap) TIOGA_FREE(sndMap); sndMap=NULL;
  if (rcvMap) TIOGA_FREE(rcvMap); rcvMap=NULL;
  //
  nsend=ns;
  nrecv=nr;
  sndMap=(int *) malloc(sizeof(int)*nsend);
  rcvMap=(int *) malloc(sizeof(int)*nrecv);
  //
  for(i=0;i<nsend;i++) sndMap[i]=snd[i];
  for(i=0;i<nrecv;i++) rcvMap[i]=rcv[i];
}

void parallelComm::getMap(int *ns, int *nr, int **snd,int **rcv)
{
  *ns=nsend;
  *nr=nrecv;
  
  *snd=sndMap;
  *rcv=rcvMap;
  return;
}
  
void parallelComm::clearPackets(PACKET *sndPack, PACKET *rcvPack)
{
  int i;
  //
  // free Send and recv data
  //
  for(i=0;i<nsend;i++)
    {
      if (sndPack[i].intData) TIOGA_FREE(sndPack[i].intData);
      if (sndPack[i].realData) TIOGA_FREE(sndPack[i].realData);
      sndPack[i].nints=sndPack[i].nreals=0;
    }
  for(i=0;i<nrecv;i++)
    {
      if (rcvPack[i].intData) TIOGA_FREE(rcvPack[i].intData);
      if (rcvPack[i].realData) TIOGA_FREE(rcvPack[i].realData);
      rcvPack[i].nints=rcvPack[i].nreals=0;
    }
  //
}


void parallelComm::initPackets(PACKET *sndPack, PACKET *rcvPack)
{
  int i;
  //
  for(i=0;i<nsend;i++)     
    {
      sndPack[i].nints=sndPack[i].nreals=0;
      sndPack[i].intData=NULL;
      sndPack[i].realData=NULL;
    }
  //
  for(i=0;i<nrecv;i++)     
    {
      rcvPack[i].nints=rcvPack[i].nreals=0;
      rcvPack[i].intData=NULL;
      rcvPack[i].realData=NULL;
    }
  //
}
