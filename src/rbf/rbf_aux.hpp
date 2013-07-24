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


#ifndef _RBF_AUX_HPP_
#define _RBF_AUX_HPP_

#include "rbf_aux.h"

namespace rbf{
   
// 
//----------------------------------------------------------
//
int rbf_minimo(int a, int b)
{
 if(a<b)
    return a;
 else
    return b; 
} 
//  
//----------------------------------------------------------
//
int rbf_maximo(int a, int b)
{
 if(a>b)
    return a;
 else
    return b; 
}
//
//----------------------------------------------------------
//
bool rbf_is_par(unsigned int x)
{
   if((x % 2)==0)
      return true;
   else
      return false;
}
//
//----------------------------------------------------------
//
bool rbf_is_impar(unsigned int x)
{
   if((x % 2)==0)
      return false;
   else
      return true;
}
//
//----------------------------------------------------------
//
template <typename T> 
T rbf_pow(T r, int beta)
{
 switch (beta) {
    case  1:
       return r;
       break;
    case  2:
       return r*r;
       break;
    case  3:
       return r*r*r;
       break;
    case  4:
       return r*r*r*r;
    case  5:
       return r*r*r*r*r;
       break;
    case  6:
       return r*r*r*r*r*r;
       break;
    case  7:
       return r*r*r*r*r*r*r;
    case  8:
       return r*r*r*r*r*r*r*r;
    default:
        return pow( r , T(beta) );
       break;
 }
}
//
//----------------------------------------------------------
//
template <typename T> 
T rbf_norm_square(const T *x, const T *y, int n)
{
   int  i;   
   T    s = 0.0;
   
   for( i=0; i < n; i++ )   
      s += ( x[i] - y[i] ) * ( x[i] - y[i] );
   
   return s;   
}   
   
} // RBF namespace

#endif //_RBF_AUX_HPP_
