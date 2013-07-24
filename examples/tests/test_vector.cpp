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

//---------------------------------------------------------
int main(int argc, char **argv)
{
  	 
  Vector<double>   A,B,C,D;
  int              n=100000;
  
  A.Reallocate(n);
  B.Reallocate(n);
  C.Reallocate(n);
  D.Reallocate(n);
  
  times.tic();
  for(int iter=0;  iter<1;  iter++)
  {
    Add(1.0,B,C);
    Add(1.0,C,D);
    /*
        cblas_daxpy(B.GetLength(),
		1.0,
		reinterpret_cast<const double*>(B.GetData()), 1,
		reinterpret_cast<double*>(C.GetData()), 1);
		
		
		        cblas_daxpy(B.GetLength(),
		1.0,
		reinterpret_cast<const double*>(B.GetData()), 1,
		reinterpret_cast<double*>(C.GetData()), 1);*/
  }  
  printf("elapsed time = %f  MFLOPS = %f\n",times.toc(),(3.0*n*n)/(times.toc()*1e9));
  


  return 0;
}
//---------------------------------------------------------





