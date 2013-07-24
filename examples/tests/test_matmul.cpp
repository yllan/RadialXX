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




int compare (const void * a, const void * b)


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
    for (i=0; i<N; i++)
    
    
    
   printf ( "  %d x %d :\n",N,N );
         
         

				//MltAdd(1.0, A, x, 0.0, y); 
		y = A*x;
		
	   
	  
             
     printf ( " \n " );          
    printf ( "  %d x %d :\n",N,N );
  


  return 0;
}
//---------------------------------------------------------




