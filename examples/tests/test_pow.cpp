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



#include "radial.h"



int compare (const void * a, const void * b)


//---------------------------------------------------------
int main(int argc, char **argv)
{
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
  
  
   ops   = 3*N;
   
   iters = 50; 
    
    
   x.Resize(N);
   y.Resize(N); 
    
//fill the mat and vec  
    
    
    
   printf ( "  %d x %d :\n",N,N );
         
         

		
		y = pow3(x);
		
	   
	  
             
     printf ( " \n " );          
    printf ( "  %d x %d :\n",N,N );
  
     for (i=0; i<5; i++)

  printf ( " \n " );
  return 0;
}
//---------------------------------------------------------




