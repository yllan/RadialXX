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


#ifndef _RBF_IMQ_DERIVATIVES_HPP
#define _RBF_IMQ_DERIVATIVES_HPP


namespace rbf{
   

//----------------------------------------------------------   
//              Data  1-D  
//----------------------------------------------------------
template <typename T>
inline T imq_1d_x(int beta , T x, T xj, T c)
{  
  T r2;
   
  r2 = (x-xj) * (x-xj);
   
  return -beta * (x-xj) * pow( r2 + c*c, -beta/2.0 - 1.0  );
}
//----------------------------------------------------------
template <typename T>
inline T imq_1d_xx(int beta , T x, T xj, T c)
{  
   T r2,c2;
   
   c2 = c*c;
   r2 = (x-xj)*(x-xj);
   
    return  -beta * pow( r2 + c2, -beta/2.0 - 1.0) * ( 1.0 + ( (-beta-2.0) * ( (x-xj)*(x-xj) ) )/(r2 + c2) );   
}
//----------------------------------------------------------   
//              Data  2-D  
//----------------------------------------------------------
 template <typename T>
inline T imq_2d_x(int beta , T x, T y, T xj, T yj, T c)
{
   T r2;
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj); 
   
   return -beta * (x-xj) * pow( r2 + c*c, -beta/2.0 - 1.0  );   
}
//----------------------------------------------------------
template <typename T>
inline T imq_2d_dy(int beta , T x, T y, T xj, T yj, T c)
{
  T r2;
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj); 
   
   return -beta * (y-yj) * pow( r2 + c*c, -beta/2.0 - 1.0  );      
   
}
//----------------------------------------------------------
template <typename T>
inline T imq_2d_xx(int beta , T x, T y, T xj, T yj, T c)
{
   T r2,c2;
   
   c2 = c*c;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj);
   
   return  -beta * pow( r2 + c2, -beta/2.0 - 1.0) * ( 1.0 + ( (-beta-2.0) * ( (x-xj)*(x-xj) ) )/(r2 + c2) );
   
}
//----------------------------------------------------------
template <typename T>
inline T imq_2d_yy(int beta , T x, T y, T xj, T yj, T c)
{
   T r2,c2;
   
   c2 = c*c;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj);
   
   return  -beta * pow( r2 + c2, -beta/2.0 - 1.0) * ( 1.0 + ( (-beta-2.0) * ( (x-xj)*(x-xj) ) )/(r2 + c2) );   
   
}
//----------------------------------------------------------
template <typename T>
inline T imq_2d_xy(int beta, T x, T y, T xj, T yj, T c)
{
   return mq_2d_xy(-1*beta , x, y, xj, yj, c);
}
//----------------------------------------------------------
template <typename T>
inline T imq_2d_yx(int beta, T x, T y, T xj, T yj, T c)
{
   return mq_2d_yz(-1*beta , x, y, xj, yj, c);
}
//----------------------------------------------------------   
//              Data  3-D  
//----------------------------------------------------------   
template <typename T> 
T imq_3d_x(int beta , T x, T y, T xj, T yj, T z, T zj, T c)
{ 
   T r2;
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj); 
   
   return -beta * (x-xj) * pow( r2 + c*c, -beta/2.0 - 1.0  );   
}
//----------------------------------------------------------
template <typename T> 
T imq_3d_y(int beta , T x, T y, T xj, T yj, T z, T zj, T c)
{ 
   T r2;
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj); 
   
   return -beta * (y-yj) * pow( r2 + c*c, -beta/2.0 - 1.0  );   
}
//----------------------------------------------------------
template <typename T> 
T imq_3d_z(int beta , T x, T y, T xj, T yj, T z, T zj, T c)
{ 
   T r2;
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj); 
   
   return -beta * (z-zj) * pow( r2 + c*c, -beta/2.0 - 1.0  );   
}
//----------------------------------------------------------
template <typename T> 
T imq_3d_xx(int beta , T x, T y, T xj, T yj, T z, T zj, T c)
{ 
   return mq_3d_xx(-1*beta , x, y, xj, yj, z, zj, c);
}
//----------------------------------------------------------
template <typename T> 
T imq_3d_yy(int beta , T x, T y, T xj, T yj, T z, T zj, T c)
{ 
   return mq_3d_yy(-1*beta , x, y, xj, yj, z, zj, c);
}
//----------------------------------------------------------
template <typename T> 
T imq_3d_zz(int beta , T x, T y, T xj, T yj, T z, T zj, T c)
{ 
   return mq_3d_zz(-1*beta , x, y, xj, yj, z, zj, c);
}

} // RBF namespace

#endif // _RBF_IMQ_DERIVATIVES_HPP_



