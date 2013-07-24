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

//#ifdef WITH_BLAS
// #define  SELDON_WITH_CBLAS
//#endif
#include "radial.h"


//----------------------------------------------------------   
template <typename T> 
T dxx(T x, T y, T xj, T yj, T c) // take 1.296757
{ 
 T r,tmp1,tmp2;   
 T beta=1.4;
   
 r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) );   
   
   if( r == 0.0)
       return 0.0;
   else {
      tmp1 =  pow(r, (T)(beta-4.0));
      tmp2 =  (x-xj) * (x-xj);
       return  ( 1.0 + beta*log(r) )  * ( pow(r, (T)(beta-2.0)) + (beta-2.0) * tmp1 * tmp2) \
      + beta * tmp1 * tmp2; 
   }   
}

//----------------------------------------------------------   
template <typename T> 
T dxx2(T x, T y, T xj, T yj, T c) // not a good idea
{ 
 T r,tmp1,tmp2;   
 T beta=1.4;
 
 if( x==xj && y==yj )
   return 0.0;
   
 r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) );   
   

  tmp1 =  pow(r, (T)(beta-4.0));
      tmp2 =  (x-xj) * (x-xj);
       return  ( 1.0 + beta*log(r) )  * ( pow(r, (T)(beta-2.0)) + (beta-2.0) * tmp1 * tmp2) \
      + beta * tmp1 * tmp2; 

}

//----------------------------------------------------------   
template <typename T> 
T dxx3(T x, T y, T xj, T yj, T c) // 
{ 
 T r,tmp1,tmp2,tmp3;   
 T beta=1.4;
 
       tmp2 =  (x-xj) * (x-xj);
    
 r = sqrt( tmp2 + (y-yj) * (y-yj) );   
   
   if( r == 0.0)
       return 0.0;
   else {
      
      tmp1 =  pow(r, beta-4.0);
      tmp3 =  pow(r, beta-2.0);
  
return ( 1.0 + beta *log(r) )  * ( tmp3 + (beta-2.0) * tmp1 * tmp2) + beta * tmp1 * tmp2; 
   }
}



//---------------------------------------------------------
void fillMatrix(Vector<double>& x, Vector<double>& y, int n)
{
   int i,j;
   
      for(i=0;i<n;i++)
      for(j=0;j<n;j++)
      {
        dxx3(x(i),y(i), x(j), y(j), 1.0);
      }
}     
//---------------------------------------------------------
int main(int argc, char **argv)
{
      timer times;
  Vector<double>   x,y;
  int              n=2000;
  int              i,j;
  int             max_iters=5;
  x.Reallocate(n);
  y.Reallocate(n);

  
  for(i=0;i<n;i++)
  {
    x(i) = rand();
    y(i) = rand();   
  }
  
  times.tic();
  for(int iter=0;  iter<max_iters;  iter++)
  { 
     fillMatrix(x,y,n);
  }  
  printf("elapsed time = %f  MFLOPS = %f\n",times.toc(), times.toc()/double(max_iters) );
  
  times.tic();
  for(int iter=0;  iter<max_iters;  iter++)
  { 
     fillMatrix(x,y,n);
  }  
  printf("elapsed time = %f  MFLOPS = %f\n",times.toc(), times.toc()/double(max_iters) );
  

  return 0;
}
//---------------------------------------------------------





