/***************************************************************************
 *
 * Copyright (C) 2009   Jose Antonio Munoz Gomez
 *        
 * This file is part of Radial++
 * http://sourceforge.net/projects/radial/
 *
 * Radial++ is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 *
 * Radial++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For more information, please see the  Home Page:
 *    http://www.dci.dgsca.unam.mx/pderbf
 *
 ***************************************************************************/


#ifdef WITH_BLAS
#define SELDON_WITH_CBLAS
#define SELDON_WITH_LAPACK
#endif

#define  SELDON_WITH_OMP
#include <omp.h>
#include "radial.h"




int compare (const void * a, const void * b){  return ( *(double*)a > *(double*)b );}


//---------------------------------------------------------
int main(int argc, char **argv)
{
  Matrix<double>      A;
  Vector<double>      x;
  Vector<double>      y;
  int                 N;
  int                i,j;
  int              iters;
  double             ops;
  double        waverage;
  double           wmaxi;
  double     wtimes[100];
  double     mflops[100];
  timer            times;
  double error; 
   N = atoi(argv[1]);
  
  
   ops   = 2.0*N*N;
   
   iters = 50; 
    
    
   A.Resize(N,N);
   x.Resize(N);
   y.Resize(N); 
    
//fill the mat and vec  
    for (i=0; i<N; i++)      for (j=0; j<N; j++)        A(i,j) = i+j;   for (i=0; i<N; i++)        x(i)= i;
    
    
    
   printf ( "  %d x %d :\n",N,N );   printf ( " \n " );   printf ( "  wtime     wmaxi     mflops \n" );
         
         
for(int itera=0; itera<10; itera++){do{	   	for(i=0;i<iters; i++)	{	   wtimes[i] = omp_get_wtime();		
				//MltAdd(1.0, A, x, 0.0, y); 
		y = A*x;
			 wtimes[i] = omp_get_wtime() - wtimes[i];	}    //sort the wtime        qsort ( (void *)wtimes, iters, sizeof(double), compare);	waverage = 0.0;	for(i= 10;i<iters-10; i++)	  waverage += wtimes[i];	  	  waverage = waverage/double(iters-10.0*2.0);	  	  wmaxi=0.0;	  for(i= 10;i<iters-10; i++)	  {	   if(wtimes[i]>wmaxi)	     wmaxi = wtimes[i];	  }
	   
	      error = (fabs(waverage-wmaxi)/fabs(waverage))*100.0;}while(  error > 5.0);     mflops[itera] = ops/(1.0E+06*waverage);     printf ( "  %.5f   %.5f    %.2f\n", waverage,wmaxi, mflops[itera]);}  printf ( " \n" );  for(int itera=0; itera<10; itera++)   printf("   mflops = %.3f\n", mflops[itera]);
             
     printf ( " \n " );          
    printf ( "  %d x %d :\n",N,N );   printf ( " \n " );
  


  return 0;
}
//---------------------------------------------------------





