//-------------------------------------------------------------------------
//
//  Copyright (C) 2009   Jose Antonio Munoz Gomez
//
//  This file is part of Radial++
//  http://sourceforge.net/projects/radial/
//
//  Radial++ is free software;  you can redistribute it and/or it under the
//  terms of the GNU Lesser General Public License as published by the Free 
//  Software Foundation; either version 3 of the License, or (at your option)
//  any later version.
//
//  Radial++ is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
//  License for more details.
//
//-------------------------------------------------------------------------


#ifndef _ODE_CONTROL_HPP_
#define _ODE_CONTROL_HPP_

namespace rbf{


//----------------------------------------------------------
//  The time stept must be satisfy
//             Tmax + stept > Tmax
//  it is similar to:
//              1   +  eps  >  1
//  so, based on Tmax we determine the h_min to satisfy Tmax+h_min>Tmax
//
template <typename T> 
T  get_hmin(T tmax)
{
  T     aux;
  T     epsi;
  bool  band=true;
  int   op,k=0;
 
   epsi = numeric_limits<T>::epsilon();

   while(band)
   {
     aux = tmax + k*epsi;
     op  = aux > tmax;
   
     if(op==1)
        band=false;
     else
        k = k +1;
    }
 
     k = k + 2;  // more restrictive

  return k*epsi;
}
//----------------------------------------------------------
template <typename T> 
bool check_stept(T tmax, T stept)
{
  T hmin = get_hmin(tmax);
  
   if(stept < hmin)
   {
      printf("!!WARNING: the initial stept = %e is too small, the simulation maybe loop infinity.\n",stept);
      printf("           the epsilon = %e\n",numeric_limits<T>::epsilon());
      printf("           the suggested stept_min = %e\n",hmin);
      return false;
   }
  return true;
}
//----------------------------------------------------------
template <typename T> 
bool correct_stept(T  t, T tmax, T &stept)
{
   T epsi = numeric_limits<T>::epsilon();
 
   if( (t + stept) > tmax )
   {
      stept = tmax - t;

      if(stept < 5.0*epsi )
      {
         printf("!!WARNING: the new stept = %e is too close to 5.0*epsi = %e\n",stept,5.0*epsi);
      }

      if(stept <= epsi)
      {
         printf("!!WARNING: the new stept = %e is less or equal to epsi = %e\n",stept,epsi);
         printf("           the loop was terminated.\n");
         printf("           the current t = %f and tmax = %e\n",t,tmax);
         return false;
      }

      if( (t + stept) <= t ) /* optional */
      {
         printf("!!ERROR: not allowed t + stept > t, stept = %e is too small.\n",stept);
         printf("         the loop was terminated.\n");
         return false;
      }
         return    true;
   }
   else{
         return true;
   }
}      

} // RBF namespace

#endif // _ODE_CONTROL_HPP_

