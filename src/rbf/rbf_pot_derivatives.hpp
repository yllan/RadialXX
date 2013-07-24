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

#ifndef _RBF_POT_DERIVATIVES_HPP_
#define _RBF_POT_DERIVATIVES_HPP_

namespace rbf{

//----------------------------------------------------------   
//              Data  1-D  
//----------------------------------------------------------
template <typename T>
inline T pot_1d_x(int beta , T x, T xj)
{  
   T r = sqrt( (x-xj)*(x-xj) );
  
   return (x-xj) * beta * pow( r , beta-2.0 ) ;
}
//----------------------------------------------------------
template <typename T>
inline T pot_1d_xx(int beta, T x, T xj)
{  
   T  r, bm2;
   
   r   = sqrt( (x-xj)*(x-xj) );
  
   bm2 = beta - 2.0;
   
   return  beta * pow( r , bm2 ) + beta * bm2 * pow( r , beta-4.0 ) * (x-xj) * (x-xj) ;  
}
//----------------------------------------------------------   
//              Data  2-D  
//----------------------------------------------------------
 template <typename T>
inline T pot_2d_x(int beta , T x, T y, T xj, T yj)
{
   T   r;
   
   r = sqrt( ( x-xj) * (x-xj) + (y-yj) * (y-yj) );
  
   return (x-xj) * beta * pow( r , beta-2.0 ) ;
}
//----------------------------------------------------------
template <typename T>
inline T pot_2d_y(int beta , T x, T y, T xj, T yj)
{
   T   r;
   
   r = sqrt( ( x-xj)*(x-xj) + (y-yj)*(y-yj) );
  
   return (y-yj) * beta * pow( r , beta-2.0 ) ;
}
//----------------------------------------------------------
template <typename T>
inline T pot_2d_xx(int beta , T x, T y, T xj, T yj)
{
   T   r, bm2;
   
   r   = sqrt( ( x-xj)*(x-xj) + (y-yj)*(y-yj) );
   
   bm2 = beta - 2.0;
   
   return  beta * pow( r , bm2 ) + beta * bm2 * pow( r , beta-4.0 ) * (x-xj) * (x-xj) ;  
}
//----------------------------------------------------------
template <typename T>
inline T pot_2d_yy(int beta , T x, T y, T xj, T yj)
{
   T   r, bm2;
    
   r = sqrt( ( x-xj)*(x-xj) + (y-yj)*(y-yj) );
   
   bm2 = beta - 2.0;

   return  beta * pow( r , bm2 ) + beta * bm2 * pow( r, beta-4.0) * (y-yj) * (y-yj) ;  
}
//----------------------------------------------------------
template <typename T>
inline T pot_2d_xy(int beta, T x, T y, T xj, T yj)
{
   T  r;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) );   

   return beta * (beta-1.0) * (x-xj) * (y-yj) * pow(r,beta-3.0) ;
}
//----------------------------------------------------------
template <typename T>
inline T pot_2d_yx(int beta, T x, T y, T xj, T yj)
{
   return pot_2d_xy(beta, x, y, xj, yj) ;
}
//----------------------------------------------------------   
//              Data  3-D  
//----------------------------------------------------------   
template <typename T> 
T pot_3d_x(int beta , T x, T y, T xj, T yj, T z, T zj)
{ 
   T   r;
   
   r = sqrt( ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj) );
  
   return (x-xj) * beta * pow( r , beta-2.0 ) ;
}
//----------------------------------------------------------
template <typename T> 
T pot_3d_y(int beta , T x, T y, T xj, T yj, T z, T zj)
{ 
   T   r;
   
   r = sqrt( ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj) );
  
   return (y-yj) * beta * pow( r , beta-2.0 ) ;
}
//----------------------------------------------------------
template <typename T> 
T pot_3d_z(int beta, T x, T y, T xj, T yj, T z, T zj)
{ 
   T   r;
   
   r = sqrt( ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj) );
  
   return (z-zj) * beta * pow( r , beta-2.0 ) ;
}
//----------------------------------------------------------
template <typename T> 
T pot_3d_xx(int beta, T x, T y, T xj, T yj, T z, T zj)
{ 
   T   r, bm2;
   
   r   = sqrt( ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj) );
   
   bm2 = beta - 2.0;
   
   return  beta * pow( r , bm2 ) + beta * bm2 * pow( r , beta-4.0 ) * (x-xj) * (x-xj) ;  
}
//----------------------------------------------------------
template <typename T> 
T pot_3d_yy(int beta, T x, T y, T xj, T yj, T z, T zj)
{ 
   T   r, bm2;
   
   r   = sqrt( ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj) );
   
   bm2 = beta - 2.0;
   
   return  beta * pow( r , bm2 ) + beta * bm2 * pow( r , beta-4.0 ) * (y-yj) * (y-yj) ;  
}
//----------------------------------------------------------
template <typename T> 
T pot_3d_zz(int beta, T x, T y, T xj, T yj, T z, T zj)
{ 
   T   r, bm2;
   
   r   = sqrt( ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj) );
   
   bm2 = beta - 2.0;
   
   return  beta * pow( r , bm2 ) + beta * bm2 * pow( r , beta-4.0 ) * (z-zj) * (z-zj) ;  
}



} // RBF namespace

#endif // _RBF_POT_DERIVATIVES_HPP_


