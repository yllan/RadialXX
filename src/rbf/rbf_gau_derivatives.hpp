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


#ifndef _RBF_GAU_DERIVATIVES_HPP_
#define _RBF_GAU_DERIVATIVES_HPP_


namespace rbf{


//----------------------------------------------------------   
//              Data  1-D  
//----------------------------------------------------------   
template <typename T>
inline T gau_1d_x(T x, T xj, T c)
{  
   T r2;
      
   r2 =   (x - xj) * (x - xj);
   
   return -2.0 * c*c * (x - xj) * exp( -r2 * c*c );
}
//----------------------------------------------------------
template <typename T>
inline T gau_1d_xx(T x, T xj, T c)
{  
   T  r2, c2;
      
   r2 =  (x - xj) * (x - xj);
   
   c2 = c * c;
   
   return 2.0 * c2 * exp( -r2 * c2 ) * ( 2.0 * c2 * (x - xj) * (x - xj) - 1.0 );
}
//----------------------------------------------------------   
//              Data  2-D  
//----------------------------------------------------------
 template <typename T>
inline T gau_2d_x(T x, T y, T xj, T yj, T c)
{
   T r2, c2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj);
   
   c2 = c * c;
   
   return -2.0 * c2 * (x - xj) * exp( -r2 * c2 );
}
//----------------------------------------------------------
template <typename T>
inline T gau_2d_y(T x, T y, T xj, T yj, T c)
{
   T r2, c2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj);
   
   c2 = c * c;
   
   return -2.0 * c2 * (y - yj) * exp( -r2 * c2 );
}
//----------------------------------------------------------
template <typename T>
inline T gau_2d_xx(T x, T y, T xj, T yj, T c)
{
   T r2,c2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj);
   
   c2 = c * c;
   
   return 2.0 * c2 * exp( -r2 * c2 ) * ( 2.0 * c2 * (x - xj) * (x - xj) - 1.0 );
}
//----------------------------------------------------------
template <typename T>
inline T gau_2d_yy(T x, T y, T xj, T yj, T c)
{
   T  r2, c2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj);
   
   c2 = c * c;
   
   return 2.0 * c2 * exp( -r2 * c2 ) * ( 2.0 * c2 * (y - yj) * (y - yj) - 1.0 ); 
}
//----------------------------------------------------------
template <typename T>
T gau_2d_xy(T x, T y, T xj, T yj, T c)
{
   T  r2, c2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj);
   
   c2 = c * c;

   return (-2.0 * c2 * (y-yj) )*(-2.0 * c2 * (x-xj) )*exp(-r2 * c2);
}
//----------------------------------------------------------
template <typename T>
T gau_2d_yx(T x, T y, T xj, T yj, T c)
{
  return gau_2d_xy(x, y, xj, yj, c);
}
//----------------------------------------------------------   
//              Data  3-D  
//----------------------------------------------------------   
template <typename T> 
T gau_3d_x(T x, T y, T xj, T yj, T z, T zj, T c)
{ 
   T r2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj);
   
   return -2.0 * c*c * (x - xj) * exp( -r2 * c*c );  
}
//----------------------------------------------------------
template <typename T> 
T gau_3d_y(T x, T y, T xj, T yj, T z, T zj, T c)
{ 
   T r2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj);
   
   return -2.0 * c*c * (y - yj) * exp( -r2 * c*c );  
}
//----------------------------------------------------------
template <typename T> 
T gau_3d_z(T x, T y, T xj, T yj, T z, T zj, T c)
{ 
   T r2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj);
   
   return -2.0 * c*c * (z - zj) * exp( -r2 * c*c );  
}
//----------------------------------------------------------
template <typename T> 
T gau_3d_xx(T x, T y, T xj, T yj, T z, T zj, T c)
{ 
   T r2,c2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj);
   
   c2 = c * c;
   
   return 2.0 * c2 * exp( -r2 * c2 ) * ( 2.0 * c2 * (x - xj) * (x - xj) - 1.0 );
}
//----------------------------------------------------------
template <typename T> 
T gau_3d_yy(T x, T y, T xj, T yj, T z, T zj, T c)
{ 
   T r2,c2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj);
   
   c2 = c * c;
   
   return 2.0 * c2 * exp( -r2 * c2 ) * ( 2.0 * c2 * (y - yj) * (y - yj) - 1.0 );
}
//----------------------------------------------------------
template <typename T> 
T gau_3d_zz(T x, T y, T xj, T yj, T z, T zj, T c)
{ 
   T r2,c2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj);
   
   c2 = c * c;
   
   return 2.0 * c2 * exp( -r2 * c2 ) * ( 2.0 * c2 * (z - zj) * (z - zj) - 1.0 );
}


} // RBF namespace

#endif // _RBF_GAU_DERIVATIVES_HPP_





