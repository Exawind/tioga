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
/**
 * Build an alternating digital tree 
 */
#include <stdio.h>
#include <stdlib.h>
#include "codetypes.h"
#include "ADT.h"
#include "buildADTrecursion.h"

void ADT::buildADT(int d, int nelements,double *elementBbox)
{
  int i,i2,j6,j,i4;
  int *elementsAvailable;
  double *adtWork;
  int adtCount,parent,level,nav;
  int side;    
  double tolerance,delta;
  FILE *fp,*fp1;
  //
  /* set dimensions and number of elements */
  //
  ndim=d;
  nelem=nelements;
  /* set element bbox pointer */
  coord=elementBbox;  
  /*
   * Allocate work arrays
   */
  elementsAvailable=(int *) malloc(sizeof(int)*nelem);
  adtWork=(double *) malloc(sizeof(double)*nelem);
  /*
   * Allocate arrays in the class
   */
  if (adtExtents) TIOGA_FREE(adtExtents);
  adtExtents=(double *) malloc(sizeof(double)*ndim);
  if (adtIntegers) TIOGA_FREE(adtIntegers);
  adtIntegers=(int *) malloc(sizeof(int)*4*nelem);
  if (adtReals) TIOGA_FREE(adtReals);
  adtReals=(double *) malloc(sizeof(double)*nelem*ndim);
  /*
   * Determine extent of elements
   */
  for(i=0;i<ndim/2;i++)
    {
      i2=2*i;
      adtExtents[i2]=BIGVALUE;
      adtExtents[i2+1]=-BIGVALUE;
   }
  for(j=0;j<nelem;j++)
   {
     j6=6*j;	 
     for(i=0;i<ndim/2;i++)
       {
	 i2=2*i;
	 adtExtents[i2]=std::min(adtExtents[i2],coord[j6+i]);
       }
       for(i=0;i<ndim/2;i++)
       {
	 i2=2*i+1;
	 adtExtents[i2]=std::max(adtExtents[i2],coord[j6+i+ndim/2]);
       }
   }
  //
  // make the extents 1% larger
  //
  tolerance=0.01;
  for(i=0;i<ndim/2;i++)
    {
      i2=2*i;
      delta=tolerance*(adtExtents[i2+1]-adtExtents[i2]);
      adtExtents[i2]-=delta;
      adtExtents[i2+1]+=delta;
    }
  //
  // Build ADT using a recursive process now
  //
  for(i=0;i<nelem;i++)
    elementsAvailable[i]=i;
  //
  // set initialvalues
  //
  adtCount=-1;
  side=0;
  parent=0;
  level=0;
  nav=nelem;
  //
  buildADTrecursion(coord,adtReals,adtWork,adtIntegers,elementsAvailable,
		    &adtCount,side,parent,level,ndim,nelem,nav);
  //TRACEI(adtCount);
  //
  // create Inverse map
  //
  //fp=fopen("adtReals.dat","w");
  //fp1=fopen("adtInts.dat","w");
  for(i=0;i<nelem;i++)
    {
      i4=4*adtIntegers[4*i];
      adtIntegers[i4+3]=i;
    }
  //for(i=0;i<nelem;i++)
  // {
  //   fprintf(fp,"%.8e %.8e %.8e %.8e %.8e %.8e\n",adtReals[6*i],adtReals[6*i+1],adtReals[6*i+2],adtReals[6*i+3],
   //                                 adtReals[6*i+4],adtReals[6*i+5]);
   // fprintf(fp1,"%d %d %d %d\n",adtIntegers[4*i],adtIntegers[4*i+1],adtIntegers[4*i+2],adtIntegers[4*i+3]);
  // }
  //fclose(fp);
  //fclose(fp1);
  TIOGA_FREE(elementsAvailable);
  TIOGA_FREE(adtWork);
}
