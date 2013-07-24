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
void test_make(int dim, int degree)
{
  Polinomio<double> pol;
   timer             time;
   int               max_iters=1;
   

   time.tic();
 for(int iter=0 ; iter < max_iters ; iter++)
 {    
          pol.make( dim, degree );
 }          
  printf("dim = %d   m = %d  time = %.3e sec\n",dim,pol.get_M(),time.toc()/double(max_iters));
  
}      
//---------------------------------------------------------
int main(int argc, char **argv)
{


  int               dim, degree;

   degree = 6;
   
      test_make(10, degree);





  return 0;
}
/*
Degree 3
dim = 1    m = 3   time = 1.010e-06 sec
dim = 2    m = 6   time = 1.810e-06 sec
dim = 3    m = 10  time = 3.030e-06 sec
dim = 4    m = 15  time = 7.950e-06 sec
dim = 5    m = 21  time = 1.589e-05 sec
dim = 6    m = 28  time = 3.632e-05 sec
dim = 7    m = 36  time = 1.087e-04 sec
dim = 8    m = 45  time = 3.456e-04 sec
dim = 9    m = 55  time = 1.205e-03 sec
dim = 10   m = 66  time = 3.856e-03 sec


Degree real = 4
dim = 1   m = 4  time = 1.300e-06 sec
dim = 2   m = 10  time = 3.200e-06 sec
dim = 3   m = 20  time = 5.800e-06 sec
dim = 4   m = 35  time = 2.070e-05 sec
dim = 5   m = 56  time = 5.250e-05 sec
dim = 6   m = 84  time = 1.977e-04 sec
dim = 7   m = 120  time = 7.482e-04 sec
dim = 8   m = 165  time = 3.315e-03 sec
dim = 9   m = 220  time = 1.552e-02 sec
dim = 10   m = 286  time = 6.579e-02 sec


degree= 5
dim = 1   m = 5  time = 1.800e-06 sec
dim = 2   m = 15  time = 5.100e-06 sec
dim = 3   m = 35  time = 1.170e-05 sec
dim = 4   m = 70  time = 3.800e-05 sec
dim = 5   m = 126  time = 1.430e-04 sec
dim = 6   m = 210  time = 6.296e-04 sec
dim = 7   m = 330  time = 3.451e-03 sec
dim = 8   m = 495  time = 1.945e-02 sec
dim = 9   m = 715  time = 1.152e-01 sec
dim = 10   m = 1001  time = 6.046e-01 sec

*/





