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


#define SELDON_WITH_CBLAS
#define SELDON_WITH_LAPACK
#include "radial.h"

/*!
 *  Tool for time measurements (borrowed from FLENS).
 */
struct timer
{
/*!
 *  Start the clock.
 */
    void tic() { _time = std::clock(); }
/*!
 *  Stop the clock.
 */ 
    double toc() const { return double(std::clock() - _time)/CLOCKS_PER_SEC; }
    std::clock_t _time;
};

timer times;

class RMatrix:public Matrix<double>
{	
 public:
   int op2;
   static int op;
   RMatrix& operator=(const RMatrix& B)
   {
   	
   }
};

int RMatrix::op = 0;
RMatrix tmp;


void add2(const RMatrix& A, const RMatrix& B, const RMatrix& C)
{
 double* a=A.GetData();
 double* b=B.GetData();
 double* c=C.GetData();
 
 for(int i=0; i<1000*1000; i++) 	
   c[i] = a[i] + b[i];
}

RMatrix& operator+(const RMatrix& A, const RMatrix& B)
{	
	printf("a+b with ops %d %d\n",A.op2,B.op2);
	
 if(RMatrix::op==0)
 {	
   tmp.Reallocate(A.GetM(),A.GetN());
   RMatrix::op=1;
   tmp.op2=tmp.op2-1;
   add2(A,B,tmp);
 } 
 else{
 	tmp.op2=tmp.op2-1;
    add2(tmp,B,tmp);
  }
  
 return tmp;
}


void add(const RMatrix& A, const RMatrix& B, const RMatrix& C)
{
 tmp.Reallocate(1000,1000);
 double* a=A.GetData();
 double* b=B.GetData();
 double* c=C.GetData();
 double* d=tmp.GetData();
 
 for(int i=0; i<1000*1000; i++) 	
   d[i] = a[i] + b[i] + c[i];
}

//---------------------------------------------------------
int main(int argc, char **argv)
{
  	 
  Matrix<double>   A,B,C,D;
  int              n=1000;
  
  A.Reallocate(n,n);
  B.Reallocate(n,n);
  C.Reallocate(n,n);
  D.Reallocate(n,n);
  
  times.tic();
  for(int iter=0;  iter<1;  iter++)
  {
    Add(1,B,C);
  }  
  printf("elapsed time = %f  MFLOPS = %f\n",times.toc(),(3.0*n*n)/(times.toc()*1e6));
  

 RMatrix  A1,B1,C1,D1,E1;  
 
  A1.Reallocate(n,n);
  B1.Reallocate(n,n);
  C1.Reallocate(n,n);
  D1.Reallocate(n,n);
  E1.Reallocate(n,n);

  B1.op2=1;
  C1.op2=2;
  D1.op2=3;
   E1.op2=4;
  times.tic();
  for(int iter=0;  iter<1;  iter++)
  {
     B1 + C1 + D1 + E1;
  }  
  printf("elapsed time = %f  MFLOPS = %f\n",times.toc(),(3.0*n*n)/(times.toc()*1e6));
  
   times.tic();
  for(int iter=0;  iter<1;  iter++)
  {
     add(B1,C1,D1);
  }  
  printf("elapsed time = %f  MFLOPS = %f\n",times.toc(),(3.0*n*n)/(times.toc()*1e6));
  

  return 0;
}
//---------------------------------------------------------





