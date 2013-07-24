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

//----------------------------------------------------------
double myf(double x, double y)
{
  return sin(x)*exp(-y) + cos(y)*exp(-x);	
}
//----------------------------------------------------------
template<typename Vec>
void make_data(int N, Vec &x, Vec &y, Vec &f)
{
  double a,b,h;
  int    cont;	

   cont = 0;
   a    = 0;
   b    = 6;
	
   x.Resize(N*N);
   y.Resize(N*N);
   f.Resize(N*N);	
	
   h = (b-a)/double(N-1);	
	
	
   for(int i=0; i<N; i++ )
   {
     for(int j=0;j<N;j++)
     {
       x(cont) = a + h*i;
       y(cont) = a + h*j;
	cont++;
     }
   }
	
	
   for(int i=0; i<N*N; i++ )
     f(i) = myf(x(i),y(i));	
	
}
