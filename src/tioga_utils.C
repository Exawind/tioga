/* This file is part of the Tioga software library */

/* Tioga  is a tool for overset grid assembly on parallel distributed systems */
/* Copyright (C) 2015 Jay Sitaraman */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA */
#include "codetypes.h"
#include "tioga_utils.h"
#include "kaiser.h"

/***
 ** find oriented bounding box for a given set of points
*/
void findOBB(double *x,double xc[3],double dxc[3],double vec[3][3],int nnodes)
{
  int i,j,k,m,i3;
  double aa[9];
  double eigenv[3];
  double trace,sume;
  double xd[3];
  double xmin[3],xmax[3];
  int nrows,ncols;
  //
  xc[0]=xc[1]=xc[2]=0;
  //
  // find centroid coordinates (not true centroid)
  //
  for(i=0;i<nnodes;i++)
    {
      i3=3*i;
      xc[0]+=x[i3];
      xc[1]+=x[i3+1];
      xc[2]+=x[i3+2];
    }
  //
  xc[0]/=nnodes;
  xc[1]/=nnodes;
  xc[2]/=nnodes;
  if (nnodes <4)
    {
      vec[0][0]=vec[1][1]=vec[2][2]=1;
      vec[0][1]=vec[1][0]=vec[1][2]=vec[2][1]=0;
      vec[0][2]=vec[2][0]=0;
      if (nnodes==1)
	{
	  dxc[0]=1e-3;
	  dxc[1]=1e-3;
	  dxc[2]=1e-3;
	  return;
	}
      else if (nnodes==2)
	{
	  dxc[0]=std::max(1e-3,fabs(x[3]-x[0]))*0.5;
	  dxc[1]=std::max(1e-3,fabs(x[4]-x[1]))*0.5;
	  dxc[2]=std::max(1e-3,fabs(x[5]-x[2]))*0.5;
          return;
	}
      else
        {
         for(i=0;i<nnodes;i++)
          {
           i3=3*i;
           for(j=0;j<3;j++)
            dxc[j]=std::max(1e-3,fabs(x[i3+j]-x[0]));
          }
	 return;
        }
    }
  //
  // find co-variance matrix
  // aa = [I11 I12 I13;I21 I22 I23;I31 I32 I33]
  //
  //
  for(i=0;i<9;i++) aa[i]=0;
  //
  for(i=0;i<nnodes;i++)
    {
      i3=3*i;
      aa[0]+=((x[i3]-xc[0])*(x[i3]-xc[0]));
      aa[4]+=((x[i3+1]-xc[1])*(x[i3+1]-xc[1]));
      aa[8]+=((x[i3+2]-xc[2])*(x[i3+2]-xc[2]));
      aa[3]+=((x[i3]-xc[0])*(x[i3+1]-xc[1]));
      aa[6]+=((x[i3]-xc[0])*(x[i3+2]-xc[2]));
      aa[7]+=((x[i3+1]-xc[1])*(x[i3+2]-xc[2]));
    }
  aa[1]=aa[3];
  aa[2]=aa[6];
  aa[5]=aa[7];
  //
  // use kaisers method to estimate
  // eigen values and vectors of the covariance matrix
  //
  nrows=3;
  ncols=3;
  kaiser(aa,nrows,ncols,eigenv,trace,sume);
  //
  // copy the eigen vector basis on to vec
  //
  m=0;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      {
       vec[i][j]=aa[m++];
      }
  //
  // find min and max bounds in the bounding box
  // vector basis
  //
  for(j=0;j<3;j++)
    {
      xmax[j]=-BIGVALUE;
      xmin[j]=BIGVALUE;
    }
  for(i=0;i<nnodes;i++)
    {
      i3=3*i;
      for(j=0;j<3;j++) xd[j]=0;
      //
      for(j=0;j<3;j++)
	for(k=0;k<3;k++)
	  xd[j]+=(x[i3+k]-xc[k])*vec[j][k];
      //
      for(j=0;j<3;j++)
	{
	  xmax[j]=std::max(xmax[j],xd[j]);
	  xmin[j]=std::min(xmin[j],xd[j]);
	}
    }
  //
  // find the extents of the box
  // and coordinates of the center w.r.t. xc
  // increase extents by 1% for tolerance
  //
  for(j=0;j<3;j++)
    {
      dxc[j]=(xmax[j]-xmin[j])*0.5*1.01;
      xd[j]=(xmax[j]+xmin[j])*0.5;
    }
  //
  // find the center of the box in
  // actual cartesian coordinates
  //
  for(j=0;j<3;j++)
    {
    for(k=0;k<3;k++)
      xc[j]+=(xd[k]*vec[k][j]);
    }
}
/**
 check if a point is inside the
 provided hole map
*/
int checkHoleMap(double *x,int *nx,int *sam,double *extents)
{
  int i;
  int mm;
  double dx[3];
  int ix[3];

  for(i=0;i<3;i++) dx[i]=(extents[i+3]-extents[i])/nx[i];
  for(i=0;i<3;i++)
    {
      ix[i]=(x[i]-extents[i])/dx[i];
      if (ix[i] < 0 || ix[i] > nx[i]-1) return 0;
    }
  mm=ix[2]*nx[1]*nx[0]+ix[1]*nx[0]+ix[0];
  return sam[mm];
}

int search_octant(double *xpt,double ds[3],ADAPTIVE_HOLEMAP *AHM,
                  octant_t *parent,int level_id){
  double dx[3];
  double xlo[3];
  uint32_t cidx;
  int e;

  const qcoord_t levelh = OCTANT_LEN(level_id);

  // compute octant lengths
  dx[0] = ds[0]*INT2DBL*levelh;
  dx[1] = ds[1]*INT2DBL*levelh;
  dx[2] = ds[2]*INT2DBL*levelh;

  // compute physical coordinates
  xlo[0] = AHM->meta.extents_lo[0] + ds[0]*INT2DBL*parent->x;
  xlo[1] = AHM->meta.extents_lo[1] + ds[1]*INT2DBL*parent->y;
  xlo[2] = AHM->meta.extents_lo[2] + ds[2]*INT2DBL*parent->z;

  /* check intersection with box around each point */
  if(xpt[0] < xlo[0])       return OUTSIDE_SB;
  if(xpt[0] > xlo[0]+dx[0]) return OUTSIDE_SB;
  if(xpt[1] < xlo[1])       return OUTSIDE_SB;
  if(xpt[1] > xlo[1]+dx[1]) return OUTSIDE_SB;
  if(xpt[2] < xlo[2])       return OUTSIDE_SB;
  if(xpt[2] > xlo[2]+dx[2]) return OUTSIDE_SB;

  // point exists inside this octant: need to check if leaf or parent
  // check if parent is refined: [yes] search children, [no] return filltype

  if(parent->leafflag == 0){
    level_t *nextLevel = &(AHM->levels[level_id+1]);

    /* Search children for filltype:
     *  Stop at first non-zero value since a point
     *  may exist in at most one leaf octant.
     *
     * Note: we set up the filltype to have any value > 0 (OUTSIDE_SB)
     *  be an INSIDE_SB point or a possible WALL_SB.
     */
    int value = OUTSIDE_SB;
    for(e=0; e<OCTANT_CHILDREN && value==OUTSIDE_SB; e++){
      cidx = parent->children[e];
      value = search_octant(xpt,ds,AHM,&(nextLevel->octants[cidx]),level_id+1);
    }
    return value;
  }

  // found point in octant which is a leaf: return fill type
  return parent->filltype;
}

int checkAdaptiveHoleMap(double *xpt,ADAPTIVE_HOLEMAP *AHM){
  double ds[3];

  // get octant physical lengths
  ds[0] = AHM->meta.extents_hi[0] - AHM->meta.extents_lo[0];
  ds[1] = AHM->meta.extents_hi[1] - AHM->meta.extents_lo[1];
  ds[2] = AHM->meta.extents_hi[2] - AHM->meta.extents_lo[2];

  // recursively search from root octant
  octant_t *root = &(AHM->levels[0].octants[0]);
  return search_octant(xpt,ds,AHM,root,0);
}

/**
 fill a given hole map using iterative
 flood fill from outside the marked boundary.
 boundary is marked by "2"
*/
void fillHoleMap(int *holeMap, int ix[3],int isym)
{
  int m;
  int ii,jj,kk,mm;
  int ns2;
  int i,j,k;
  int ipaint,npaint,nneig;
  int mk[6];
  //
  // now start from outside and paint the
  // the exterior
  //
  ns2=ix[0]*ix[1];
  //
  for(kk=0;kk<ix[2];kk+=(ix[2]-1))
    for(jj=0;jj<ix[1];jj++)
      for(ii=0;ii<ix[0];ii++)
        {
	  mm=kk*ns2+jj*ix[0]+ii;
	  holeMap[mm]=1;
	}
  for(kk=0;kk<ix[2];kk++)
    for(jj=0;jj<ix[1];jj+=(ix[1]-1))
      for(ii=0;ii<ix[0];ii++)
        {
  	  mm=kk*ns2+jj*ix[0]+ii;
         holeMap[mm]=1;
  	}
  for(kk=0;kk<ix[2];kk++)
    for(jj=0;jj<ix[1];jj++)
      for(ii=0;ii<ix[0];ii+=(ix[0]-1))
        {
	  mm=kk*ns2+jj*ix[0]+ii;
	  holeMap[mm]=1;
	}
  npaint=ns2*ix[2];
  while(npaint > 0)
    {
      npaint=0;
      for(k=1;k<ix[2]-1;k++)
        for(j=1;j<ix[1]-1;j++)
          for(i=1;i<ix[0]-1;i++)
            {
              m=k*ns2+j*ix[0]+i;
              if (holeMap[m]==0)
                {
                  ipaint=0;
                  if (isym==1)
                   {
		     mk[0]=m-ns2;
		     mk[1]=m+ns2;
		     mk[2]=m-ix[0];
		     mk[3]=m+ix[0];
		     mk[4]=m-1;
		     mk[5]=m+1;
		     nneig=4;
		   }
		  else if (isym==2)
		    {
		     mk[0]=m-ns2;
		     mk[1]=m+ns2;
		     mk[4]=m-ix[0];
		     mk[5]=m+ix[0];
		     mk[2]=m-1;
		     mk[3]=m+1;
		     nneig=4;
		    }
		  else if (isym==3)
		    {
		     mk[4]=m-ns2;
		     mk[5]=m+ns2;
		     mk[0]=m-ix[0];
		     mk[1]=m+ix[0];
		     mk[2]=m-1;
		     mk[3]=m+1;
		     nneig=4;
		    }
		  else
		    {
		      mk[0]=m-ns2;
		      mk[1]=m+ns2;
		      mk[2]=m-ix[0];
		      mk[3]=m+ix[0];
		      mk[4]=m-1;
		      mk[5]=m+1;
		      nneig=6;
		    }
                  for (kk=0;kk<nneig && ipaint==0;kk++)
		    {
		      ipaint=(ipaint || holeMap[mk[kk]]==1);
		    }
                  if (ipaint > 0)
                    {
                      holeMap[m]=1;
                      npaint++;
                    }
                }
            }
    }
  for(i=0;i<ix[2]*ix[1]*ix[0];i++)
   {
    holeMap[i]=(holeMap[i] ==0 || holeMap[i]==2);
   }
}

int obbIntersectCheck(double vA[3][3],double xA[3],double dxA[3],
		       double vB[3][3],double xB[3],double dxB[3])
{
  int iflag;
  int i,j,k;
  int i1,i2,j1,j2;
  double r,r0,r1;
  double d1,d2;
  double eps=1e-12;
  double D[3];
  double c[3][3];
  //
  // D=distance between centers
  // C=scalar product of axes
  //
  for(i=0;i<3;i++) D[i]=xB[i]-xA[i];
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      {
	c[i][j]=0;
	for(k=0;k<3;k++)
	  c[i][j]=c[i][j]+vA[i][k]*vB[j][k];
      }
  //
  // separating axes based on the faces of box A
  //
  for(i=0;i<3;i++)
    {
      r0=dxA[i];
      r1=0;
      r=0;
      for(j=0;j<3;j++)
	{
	  r1+=dxB[j]*fabs(c[i][j]);
	  r+=fabs(vA[i][j])*D[j];
	}
      if (r > (r0+r1+eps)) return 0;
    }
  //
  // separating axes based on the faces of box B
  //
  for(i=0;i<3;i++)
    {
      r1=dxB[i];
      r0=0;
      r=0;
      for(j=0;j<3;j++)
	{
	  r0+=dxA[j]*fabs(c[j][i]);
	  r+=fabs(vB[i][j])*D[j];
	}
      if (r > (r0+r1+eps)) return 0;
    }
  //
  // cross products
  //
  for(i=0;i<3;i++)
    {
      i1=(i+1)%3;
      i2=(i+2)%3;
      for(j=0;j<3;j++)
	{
	  j1=(j+1)%3;
	  j2=(j+2)%3;

	  r0=dxA[i1]*fabs(c[i2][j])+dxA[i2]*fabs(c[i1][j]);
	  r1=dxB[j1]*fabs(c[i][j2])+dxB[j2]*fabs(c[i][j1]);

	  d2=0;
	  d1=0;
	  for(k=0;k<3;k++)
	    {
	      d2+=vA[i2][k]*D[k];
	      d1+=vA[i1][k]*D[k];
	    }

	  r=fabs(c[i1][j]*d2-c[i2][j]*d1);

	  if (r > (r0+r1+eps)) {
	    return 0;
	  }
	}
    }
  //
  // return zero if no separation can be found
  //
  return 1;
}

void getobbcoords(double xc[3],double dxc[3],double vec[3][3],double xv[8][3])
{
  int i,j,k,ik;
  for(i=0;i<8;i++)
   {
     for(k=0;k<3;k++)
      xv[i][k]=xc[k];
     for(k=0;k<3;k++)
       {
         ik=(2*((i & (1 << k)) >> k)-1);
         for(j=0;j<3;j++)
           xv[i][j]+=ik*dxc[k]*vec[k][j];
       }
   }
}

void transform2OBB(double xv[3],double xc[3],double vec[3][3],double xd[3])
{
 int j,k;
 for(j=0;j<3;j++)
  {
   xd[j]=0;
   for(k=0;k<3;k++)
      xd[j]+=(xv[k]-xc[k])*vec[j][k];
  }
}


void writebbox(OBB *obb,int bid)
{
  FILE *fp;
  char intstring[12];
  char fname[80];
  int l,k,j,m,il,ik,ij;
  REAL xx[3];

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"qbox%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Box file\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEPOINT\n",8,
	  1);

  for(l=0;l<2;l++)
    {
      il=2*(l%2)-1;
      for(k=0;k<2;k++)
	{
	  ik=2*(k%2)-1;
	  for(j=0;j<2;j++)
	    {
	      ij=2*(j%2)-1;
	      xx[0]=xx[1]=xx[2]=0;
	      for(m=0;m<3;m++)
		xx[m]=obb->xc[m]+ij*obb->vec[0][m]*obb->dxc[0]
		  +ik*obb->vec[1][m]*obb->dxc[1]
		  +il*obb->vec[2][m]*obb->dxc[2];
	      fprintf(fp,"%f %f %f\n",xx[0],xx[1],xx[2]);
	    }
	}
    }
  fprintf(fp,"1 2 4 3 5 6 8 7\n");
  fclose(fp);
}
void writebboxdiv(OBB *obb,int bid)
{
  double mapdx[3],mdx[3],mx0[3];
  double mapdims[3];
  int ncells,npts;
  double xv[8][3],xc[3],xd[3];
  int i,j,k,l,m,n;
  int iorder[8]={1, 2, 4, 3, 5, 6, 8, 7};
  FILE *fp;
  char intstring[12];
  char fname[80];

  for(j=0;j<3;j++) { mapdims[j]=12; mapdx[j]=2*obb->dxc[j]/mapdims[j]; mdx[j]=0.5*mapdx[j];mx0[j]=0;}
  ncells=mapdims[2]*mapdims[1]*mapdims[0];
  npts=ncells*8;
  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"dbox%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Box file\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\"\n");
  fprintf(fp,"ZONE T=\"VOL_MIXED\",N=%d E=%d ET=BRICK, F=FEPOINT\n",npts,ncells);

  getobbcoords(mx0,mdx,obb->vec,xv);
  for(l=0;l<mapdims[2];l++)
    for(k=0;k<mapdims[1];k++)
      for(j=0;j<mapdims[0];j++)
        {
          xd[0]=-obb->dxc[0]+j*mapdx[0]+mapdx[0]*0.5;
          xd[1]=-obb->dxc[1]+k*mapdx[1]+mapdx[1]*0.5;
          xd[2]=-obb->dxc[2]+l*mapdx[2]+mapdx[2]*0.5;
          for(n=0;n<3;n++)
            {
              xc[n]=obb->xc[n];
              for(m=0;m<3;m++)
                xc[n]+=(xd[m]*obb->vec[m][n]);
            }
          for(m=0;m<8;m++)
            fprintf(fp,"%f %f %f\n",xc[0]+xv[m][0],xc[1]+xv[m][1],xc[2]+xv[m][2]);
        }
  for(l=0;l<ncells;l++) {
    for(m=0;m<8;m++)
      fprintf(fp,"%d ",iorder[m]+l*8);
    fprintf(fp,"\n");
  }
}


void writePoints(double *x,int nsearch,int bid)
{
  FILE *fp;
  char intstring[12];
  char fname[80];
  int i;

  sprintf(intstring,"%d",100000+bid);
  sprintf(fname,"points%s.dat",&(intstring[1]));
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"Box file\"\n");
  fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\"\n");
  for(i=0;i<nsearch;i++)
    fprintf(fp,"%f %f %f\n",x[3*i],x[3*i+1],x[3*i+2]);
  fclose(fp);
}
/*
 * Create a unique hash for list of coordinates with duplicates in
 * them. Find the rtag as max of all duplicate samples. itag contains
 * the hash to the real point
 */
void uniquenodes(double *x,int *meshtag,double *rtag,int *itag,int *nn)
{
  int NSUB=101;
  int i,j,k,m,ij,i3,jj,kk,ll,p1,p2,indx,jmax,kmax,lmax,nsblks,jkmax;
  double xmax[3],xmin[3],ds,dsi,dsx,dsxi,dsy,dsyi,dsz,dszi;
  int *cft,*numpts,*ilist;
  int nnodes=*nn;

  for(j=0;j<3;j++) xmax[j]=-1E15;
  for(j=0;j<3;j++) xmin[j]=1E15;

  for(i=0;i<nnodes;i++)
    for(j=0;j<3;j++) {
      xmax[j]=std::max(xmax[j],x[3*i+j]);
      xmin[j]=std::min(xmin[j],x[3*i+j]);
    }

  ds=(xmax[0]-xmin[0]+xmax[1]-xmin[1]+xmax[2]-xmin[2])/3.0/NSUB;
  dsi=1.0/ds;
  for(j=0;j<3;j++) xmax[j]+=ds;
  for(j=0;j<3;j++) xmin[j]-=ds;
  jmax=std::min(round((xmax[0]-xmin[0])*dsi),static_cast<double>(NSUB));
  jmax=std::max(jmax,1);
  dsx=(xmax[0]-xmin[0]+TOL)/jmax;
  dsxi=1./dsx;
  kmax=std::min(round((xmax[1]-xmin[1])*dsi),static_cast<double>(NSUB));
  kmax=std::max(kmax,1);
  dsy=(xmax[1]-xmin[1]+TOL)/kmax;
  dsyi=1./dsy;
  lmax=std::min(round((xmax[2]-xmin[2])*dsi),static_cast<double>(NSUB));
  lmax=std::max(lmax,1);
  dsz=(xmax[2]-xmin[2]+TOL)/lmax;
  dszi=1./dsz;
  nsblks=jmax*kmax*lmax;
  jkmax=jmax*kmax;
  cft=(int *)malloc(sizeof(int)*(nsblks+1));
  numpts=(int *)malloc(sizeof(int)*nsblks);
  ilist=(int *)malloc(sizeof(int)*nnodes);

  for(i=0;i<nsblks;i++) numpts[i]=0;
  for(i=0;i<nnodes;i++)
    {
      i3=3*i;
      jj=(int)((x[i3]-xmin[0])*dsxi);
      kk=(int)((x[i3+1]-xmin[1])*dsyi);
      ll=(int)((x[i3+2]-xmin[2])*dszi);
      indx=ll*jkmax+kk*jmax+jj;
      numpts[indx]=numpts[indx]+1;
    }

  cft[0]=0;
  for(i=0;i<nsblks;i++) cft[i+1]=cft[i]+numpts[i];

  for(i=0;i<nnodes;i++)
    {
      i3=3*i;
      jj=(int)((x[i3]-xmin[0])*dsxi);
      kk=(int)((x[i3+1]-xmin[1])*dsyi);
      ll=(int)((x[i3+2]-xmin[2])*dszi);
      indx=ll*jkmax+kk*jmax+jj;
      ilist[cft[indx]+numpts[indx]-1]=i;
      numpts[indx]--;
      itag[i]=i;
    }

  for(i=0;i<nsblks;i++)
    for(j=cft[i];j<cft[i+1];j++)
      {
	p1=ilist[j];
	for(k=j+1;k<cft[i+1];k++)
	  {
	    p2=ilist[k];
	    if ( fabs(x[3*p1  ]-x[3*p2  ])+
		 fabs(x[3*p1+1]-x[3*p2+1])+
		 fabs(x[3*p1+2]-x[3*p2+2]) < TOL &&
                 meshtag[p1]==meshtag[p2])
	      {
		if (p1 > p2) {
		  rtag[p2]=std::max(rtag[p1],rtag[p2]);
		  itag[p1]=itag[p2];
		}
		else {
		  rtag[p1]=std::max(rtag[p1],rtag[p2]);
		  itag[p2]=itag[p1];
		}
	      }
	  }
      }
  /*
  m=0;
  for(i=0;i<nnodes;i++)
    if (itag[i]==i+1) {
     m++;
   }
  */
  TIOGA_FREE(ilist);
  TIOGA_FREE(cft);
  TIOGA_FREE(numpts);
}
//
// modify ADT builder to remove common nodes
//
void uniqNodesTree(double *coord,
		   int *itag,double *rtag,int *meshtag,
		   int *elementsAvailable,
		   int ndim, int nav)
{
  int nd=ndim;
  double coordmid;
  int i,j,ibox;
  int p1,p2;
  int *tmpint;
  int npts[8],iv[3],cft[9];
  double xmax[3],xmin[3],xmid[3],dx[3],xp[3];
  int icheck=1;
  //
  // if there are more than 10 elements divide the tree
  //
  if (nav > 20) {
    //
    // find the bound of the boxes
    //
    icheck=0;
    xmin[0]=xmin[1]=xmin[2]=BIGVALUE;
    xmax[0]=xmax[1]=xmax[2]=-BIGVALUE;
    for(i=0;i<nav;i++)
      for(j=0;j<nd;j++)
	{
	  xmin[j]=std::min(xmin[j],coord[ndim*elementsAvailable[i]+j]);
	  xmax[j]=std::max(xmax[j],coord[ndim*elementsAvailable[i]+j]);
	}
    for(j=0;j<nd;j++) {
      xmid[j]=(xmax[j]+xmin[j])*0.5;
      dx[j]=(xmax[j]-xmin[j])*0.5+TOL;
    }
    for(j=0;j<8;j++) npts[j]=0;
    for(i=0;i<nav;i++)
      {
	for(j=0;j<3;j++) {
	  xp[j]=coord[ndim*elementsAvailable[i]+j]-xmid[j];
	  iv[j]=floor(xp[j]/dx[j])+1;
	}
	ibox= 4*iv[0]+2*iv[1]+iv[2];
	npts[ibox]++;
      }
    for(j=0;j<8;j++) if(npts[j]==nav) icheck=1;
    if (!icheck) {
    cft[0]=0;
    for(j=0;j<8;j++)
      cft[j+1]=cft[j]+npts[j];
    tmpint=(int *)malloc(sizeof(int)*nav);
    for(i=0;i<nav;i++)
      {
	for(j=0;j<3;j++){
	  xp[j]=coord[ndim*elementsAvailable[i]+j]-xmid[j];
	  iv[j]=floor(xp[j]/dx[j])+1;
	}
	ibox= 4*iv[0]+2*iv[1]+iv[2];
	npts[ibox]=npts[ibox]-1;
	tmpint[npts[ibox]+cft[ibox]]=elementsAvailable[i];
      }
    for(i=0;i<nav;i++)
      elementsAvailable[i]=tmpint[i];
    TIOGA_FREE(tmpint);
    for(j=0;j<8;j++)
      if (cft[j+1] > cft[j])
        uniqNodesTree(coord,itag,rtag,meshtag,&(elementsAvailable[cft[j]]),
	    ndim,cft[j+1]-cft[j]);
      }
  }
  if (icheck) {
    for(i=0;i<nav;i++)
      {
	p1=elementsAvailable[i];
	for(j=i+1;j<nav;j++)
	  {
	    p2=elementsAvailable[j];
	    if (fabs(coord[3*p1  ]-coord[3*p2  ])+
		fabs(coord[3*p1+1]-coord[3*p2+1])+
		fabs(coord[3*p1+2]-coord[3*p2+2]) < TOL &&
		meshtag[p1]==meshtag[p2])
	      {
		if (p1 > p2) {
		  rtag[p2]=std::max(rtag[p1],rtag[p2]);
		  itag[p1]=itag[p2];
		}
		else {
		  rtag[p1]=std::max(rtag[p1],rtag[p2]);
		  itag[p2]=itag[p1];
		}
	      }
	  }
      }
  }
}
/*
 * Create a unique hash for list of coordinates with duplicates in
 * them. Find the rtag as max of all duplicate samples. itag contains
 * the hash to the real point
 */
void uniquenodes_octree(double *x,int *meshtag,double *rtag,int *itag,
		     int *nn)
{
  int nelem=*nn;
  int *elementsAvailable;
  int i;

  elementsAvailable=(int *)malloc(sizeof(int)*nelem);
  for(i=0;i<nelem;i++)
    {
     elementsAvailable[i]=i;
     itag[i]=i;
    }

  uniqNodesTree(x,itag,rtag,meshtag,elementsAvailable,3,nelem);

  TIOGA_FREE(elementsAvailable);
}

void qcoord_to_vertex(qcoord_t x,
                      qcoord_t y,
                      qcoord_t z,
                      double *vertices,
                      double vxyz[3]){

  int xi, yi, zi;
  double wx[2], wy[2], wz[2];
  double xfactor, yfactor;

  vxyz[0] = vxyz[1] = vxyz[2] = 0.;

  wx[1] = (double) x / (double) OCTANT_ROOT_LEN;
  wx[0] = 1. - wx[1];

  wy[1] = (double) y / (double) OCTANT_ROOT_LEN;
  wy[0] = 1. - wy[1];

  wz[1] = (double) z / (double) OCTANT_ROOT_LEN;
  wz[0] = 1. - wz[1];

  int vindex = 0;
  for (zi = 0; zi < 2; ++zi) {
    for (yi = 0; yi < 2; ++yi) {
      yfactor = wz[zi] * wy[yi];
      for (xi = 0; xi < 2; ++xi) {
        xfactor = yfactor * wx[xi];

        vxyz[0] += xfactor * vertices[3*vindex + 0];
        vxyz[1] += xfactor * vertices[3*vindex + 1];
        vxyz[2] += xfactor * vertices[3*vindex + 2];
        vindex++;
      }
    }
  }
}

uint8_t octant_filltype(octant_full_t *c, const qcoord_t inc){
  /* flood-fill octant boundaries to (OUTSIDE_SB) if touching
   * the global boundary; used to initialize the flood fill process.
   */
  // returns (OUTSIDE_SB) if touching an octant boundary, (INSIDE_SB) otherwise
  uint8_t type = (((c->x == 0) || (c->y == 0) || (c->z == 0) ||
                  ((c->x + inc) == OCTANT_ROOT_LEN) ||
                  ((c->y + inc) == OCTANT_ROOT_LEN) ||
                  ((c->z + inc) == OCTANT_ROOT_LEN))) ? OUTSIDE_SB:INSIDE_SB;
  return type;
}

void octant_children(uint8_t children_level,uint32_t newidx,
                     octant_full_t * q,
                     octant_full_t *c0, octant_full_t *c1,
                     octant_full_t *c2, octant_full_t *c3,
                     octant_full_t *c4, octant_full_t *c5,
                     octant_full_t *c6, octant_full_t *c7){
  /* FILL IN OCTANT CHILDREN DATA:
   * 0. valid flag
   * 1. coordinates
   * 2. neighboring octants
   * 3. initial fill type
   */

  const qcoord_t inc = OCTANT_LEN (children_level);

  // set parent's children octant pointers
  q->children[0] = c0;
  q->children[1] = c1;
  q->children[2] = c2;
  q->children[3] = c3;
  q->children[4] = c4;
  q->children[5] = c5;
  q->children[6] = c6;
  q->children[7] = c7;

  // set children information
  // LO-LO-LO OCTANT
  c0->x = q->x;
  c0->y = q->y;
  c0->z = q->z;
  c0->id = newidx+0;
  c0->filltype = octant_filltype(c0,inc);
  c0->refined = 0;

  // HI-LO-LO OCTANT
  c1->x = c0->x | inc;
  c1->y = c0->y;
  c1->z = c0->z;
  c1->id = newidx+1;
  c1->filltype = octant_filltype(c1,inc);
  c1->refined = 0;

  // LO-HI-LO OCTANT
  c2->x = c0->x;
  c2->y = c0->y | inc;
  c2->z = c0->z;
  c2->id = newidx+2;
  c2->filltype = octant_filltype(c2,inc);
  c2->refined = 0;

  // HI-HI-LO OCTANT
  c3->x = c1->x;
  c3->y = c2->y;
  c3->z = c0->z;
  c3->id = newidx+3;
  c3->filltype = octant_filltype(c3,inc);
  c3->refined = 0;

  // LO-LO-HI OCTANT
  c4->x = c0->x;
  c4->y = c0->y;
  c4->z = c0->z | inc;
  c4->id = newidx+4;
  c4->filltype = octant_filltype(c4,inc);
  c4->refined = 0;

  // HI-LO-HI OCTANT
  c5->x = c1->x;
  c5->y = c1->y;
  c5->z = c4->z;
  c5->id = newidx+5;
  c5->filltype = octant_filltype(c5,inc);
  c5->refined = 0;

  // LO-HI-HI OCTANT
  c6->x = c2->x;
  c6->y = c2->y;
  c6->z = c4->z;
  c6->id = newidx+6;
  c6->filltype = octant_filltype(c6,inc);
  c6->refined = 0;

  // HI-HI-HI OCTANT
  c7->x = c3->x;
  c7->y = c3->y;
  c7->z = c4->z;
  c7->id = newidx+7;
  c7->filltype = octant_filltype(c7,inc);
  c7->refined = 0;
}

void octant_children_neighbors(const octant_full_t * q,
                               octant_full_t *c0, octant_full_t *c1,
                               octant_full_t *c2, octant_full_t *c3,
                               octant_full_t *c4, octant_full_t *c5,
                               octant_full_t *c6, octant_full_t *c7){
  // 3D OCTANTS: lexicographic order [x][y][z]
  // (Drawn as 2D: 0,1,2,3 octants are zlo & 4,5,6,7 octants are zhi)
  /*  *-------*-------*  *-------*-------*       y
   *  |       |       |  |       |       |       ^
   *  |   2   |   3   |  |   6   |   7   |       |
   *  |       |       |  |       |       |       |
   *  *-------*-------*  *-------*-------*       *-----> x
   *  |       |       |  |       |       |        \
   *  |   0   |   1   |  |   4   |   5   |         \
   *  |       |       |  |       |       |          \
   *  *-------*-------*  *-------*-------*           z
   */

  /* This function assigns each child octant its neighbors.
   * A neighbor may come from outside the original parent footprint,
   * thus we use the parent's neighbor to determine the child's neighbor.
   * In this event, we must check if the parent's neighbor exists,
   * e.g. q->nhbr[XLO] might be NULL, and if it exists, check if it is refined.
   * If so, the neighbor on this level will be a child of the parent's neighbor.
   * Note we take advantage of short circuiting in the ternary logic check.
   */

  // OCTANT: [0]
  c0->nhbr[XLO] = (q->nhbr[XLO] && q->nhbr[XLO]->refined) ? q->nhbr[XLO]->children[1]:q->nhbr[XLO];
  c0->nhbr[XHI] = c1;
  c0->nhbr[YLO] = (q->nhbr[YLO] && q->nhbr[YLO]->refined) ? q->nhbr[YLO]->children[2]:q->nhbr[YLO];
  c0->nhbr[YHI] = c2;
  c0->nhbr[ZLO] = (q->nhbr[ZLO] && q->nhbr[ZLO]->refined) ? q->nhbr[ZLO]->children[4]:q->nhbr[ZLO];
  c0->nhbr[ZHI] = c4;

  // OCTANT: [1]
  c1->nhbr[XLO] = c0;
  c1->nhbr[XHI] = (q->nhbr[XHI] && q->nhbr[XHI]->refined) ? q->nhbr[XHI]->children[0]:q->nhbr[XHI];
  c1->nhbr[YLO] = (q->nhbr[YLO] && q->nhbr[YLO]->refined) ? q->nhbr[YLO]->children[3]:q->nhbr[YLO];
  c1->nhbr[YHI] = c3;
  c1->nhbr[ZLO] = (q->nhbr[ZLO] && q->nhbr[ZLO]->refined) ? q->nhbr[ZLO]->children[5]:q->nhbr[ZLO];
  c1->nhbr[ZHI] = c5;

  // OCTANT: [2]
  c2->nhbr[XLO] = (q->nhbr[XLO] && q->nhbr[XLO]->refined) ? q->nhbr[XLO]->children[3]:q->nhbr[XLO];
  c2->nhbr[XHI] = c3;
  c2->nhbr[YLO] = c0;
  c2->nhbr[YHI] = (q->nhbr[YHI] && q->nhbr[YHI]->refined) ? q->nhbr[YHI]->children[0]:q->nhbr[YHI];
  c2->nhbr[ZLO] = (q->nhbr[ZLO] && q->nhbr[ZLO]->refined) ? q->nhbr[ZLO]->children[6]:q->nhbr[ZLO];
  c2->nhbr[ZHI] = c6;

  // OCTANT: [3]
  c3->nhbr[XLO] = c2;
  c3->nhbr[XHI] = (q->nhbr[XHI] && q->nhbr[XHI]->refined) ? q->nhbr[XHI]->children[2]:q->nhbr[XHI];
  c3->nhbr[YLO] = c1;
  c3->nhbr[YHI] = (q->nhbr[YHI] && q->nhbr[YHI]->refined) ? q->nhbr[YHI]->children[1]:q->nhbr[YHI];
  c3->nhbr[ZLO] = (q->nhbr[ZLO] && q->nhbr[ZLO]->refined) ? q->nhbr[ZLO]->children[7]:q->nhbr[ZLO];
  c3->nhbr[ZHI] = c7;

  // OCTANT: [4]
  c4->nhbr[XLO] = (q->nhbr[XLO] && q->nhbr[XLO]->refined) ? q->nhbr[XLO]->children[5]:q->nhbr[XLO];
  c4->nhbr[XHI] = c5;
  c4->nhbr[YLO] = (q->nhbr[YLO] && q->nhbr[YLO]->refined) ? q->nhbr[YLO]->children[6]:q->nhbr[YLO];
  c4->nhbr[YHI] = c6;
  c4->nhbr[ZLO] = c0;
  c4->nhbr[ZHI] = (q->nhbr[ZHI] && q->nhbr[ZHI]->refined) ? q->nhbr[ZHI]->children[0]:q->nhbr[ZHI];

  // OCTANT: [5]
  c5->nhbr[XLO] = c4;
  c5->nhbr[XHI] = (q->nhbr[XHI] && q->nhbr[XHI]->refined) ? q->nhbr[XHI]->children[4]:q->nhbr[XHI];
  c5->nhbr[YLO] = (q->nhbr[YLO] && q->nhbr[YLO]->refined) ? q->nhbr[YLO]->children[7]:q->nhbr[YLO];
  c5->nhbr[YHI] = c7;
  c5->nhbr[ZLO] = c1;
  c5->nhbr[ZHI] = (q->nhbr[ZHI] && q->nhbr[ZHI]->refined) ? q->nhbr[ZHI]->children[1]:q->nhbr[ZHI];

  // OCTANT: [6]
  c6->nhbr[XLO] = (q->nhbr[XLO] && q->nhbr[XLO]->refined) ? q->nhbr[XLO]->children[7]:q->nhbr[XLO];
  c6->nhbr[XHI] = c7;
  c6->nhbr[YLO] = c4;
  c6->nhbr[YHI] = (q->nhbr[YHI] && q->nhbr[YHI]->refined) ? q->nhbr[YHI]->children[4]:q->nhbr[YHI];
  c6->nhbr[ZLO] = c2;
  c6->nhbr[ZHI] = (q->nhbr[ZHI] && q->nhbr[ZHI]->refined) ? q->nhbr[ZHI]->children[2]:q->nhbr[ZHI];

  // OCTANT: [7]
  c7->nhbr[XLO] = c6;
  c7->nhbr[XHI] = (q->nhbr[XHI] && q->nhbr[XHI]->refined) ? q->nhbr[XHI]->children[6]:q->nhbr[XHI];
  c7->nhbr[YLO] = c5;
  c7->nhbr[YHI] = (q->nhbr[YHI] && q->nhbr[YHI]->refined) ? q->nhbr[YHI]->children[5]:q->nhbr[YHI];
  c7->nhbr[ZLO] = c3;
  c7->nhbr[ZHI] = (q->nhbr[ZHI] && q->nhbr[ZHI]->refined) ? q->nhbr[ZHI]->children[3]:q->nhbr[ZHI];
}

void floodfill_octant(octant_full_t *o){
  const int nneig = 6;
  char ipaint;
  int n;

  // return if outside or wall SB
  if(o->filltype==OUTSIDE_SB ||
     o->filltype==WALL_SB) return;

  // INSIDE_SB: search neighbors
  ipaint = INSIDE_SB;
  for(n=0; n<nneig && ipaint==INSIDE_SB; n++){
    // search neighbor if it exists and is a leaf
    if(o->nhbr[n] && !o->nhbr[n]->refined){
      if(o->nhbr[n]->filltype==OUTSIDE_SB) ipaint = OUTSIDE_SB;
    }
  }
  o->filltype = ipaint;
}

/* recursive search in all directions */
void floodfill_level(level_octant_t *level){
  int nneig = 6;
  int j,n;

  for(j=0;j<level->elem_count;j++){
    octant_full_t *o = &level->octants[j];
    if(!o->refined) floodfill_octant(o);

    // fill neighbors (required since we jump around mesh due to Morton order)
    for(n=0; n<nneig; n++){
      if(o->nhbr[n] && !o->nhbr[n]->refined) floodfill_octant(o->nhbr[n]);
    }
  }
}

char checkFaceBoundaryNodes(int *nodes,const char *bcnodeflag,
                            const int numfaceverts,const int *faceConn,
                            const char *duplicatenodeflag){
  /* NOTE: faceConn is base 1 but nodes is base 0 */
  char bcFlag = 1;
  int v;

  // if any node on face is duplicate tag (wall + outer bc), return 0 (false)
  if(duplicatenodeflag){
    for(v=0;v<numfaceverts;v++){
      if(duplicatenodeflag[nodes[faceConn[v]-BASE]]){
        printf("[tioga] WARNING: Duplicate tag found -- disabling overset node %d\n",
                nodes[faceConn[v]-BASE]);
        return 0;
      }
    }
  }

  for(v=0;v<numfaceverts;v++) bcFlag &= bcnodeflag[nodes[faceConn[v]-BASE]];
  return bcFlag;
}