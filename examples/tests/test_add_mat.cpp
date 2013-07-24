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


enum ExprType {a,b,c};

#include "radial.h"

 #include <sys/time.h>
 
struct timeval posix_time1;
struct timeval posix_time2;

void tic()
{
  gettimeofday(&posix_time1, 0);
}  
double toc(void)
{
        gettimeofday(&posix_time2, 0);
        
        const double tmp_time1 = posix_time1.tv_sec + posix_time1.tv_usec * 1.0e-6;
        const double tmp_time2 = posix_time2.tv_sec + posix_time2.tv_usec * 1.0e-6;
        
        return tmp_time2 - tmp_time1;
}        
      
//---------------------------------------------------------
int main(int argc, char **argv)
{
  Matrix<double>  A,B,C,D,Z;
  timer          time;
  int             n = 100;
  int            max_iters=200;
  
  A.Resize(n,n);
  B.Resize(n,n);
  C.Resize(n,n);
  D.Resize(n,n);
  Z.Resize(n,n);
  
  A.Fill(1);
  B.Fill(2);
  C.Fill(3);
  D.Fill(4);
  
 for(int iter=0 ; iter < max_iters ; iter++)
 {    
           Z = A + B + C + D;
 }     
     
   tic();
 for(int iter=0 ; iter < max_iters ; iter++)
 {    
           Z = A + B + C + D;
 }          
  printf("Elapsed time = %e sec\n",toc()/double(max_iters));
  
     
   time.tic();
 for(int iter=0 ; iter < max_iters ; iter++)
 {    
           Z = A + B + C + D;
 }          
  printf("Elapsed time = %e sec\n",time.toc()/double(max_iters));
  



  return 0;
}
//---------------------------------------------------------





