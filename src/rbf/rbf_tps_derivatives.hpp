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


#ifndef _RBF_TPS_DERIVATIVES_HPP_
#define _RBF_TPS_DERIVATIVES_HPP_



namespace rbf{
   

//----------------------------------------------------------   
//              Data  1-D  
//----------------------------------------------------------   
template <typename T> 
T tps_1d_x(int beta, T x, T xj)
{ 
   T r;   
   
   r = sqrt( (x-xj) * (x-xj) );   
   
   if( r == 0)
       return 0.0;
   else 
       return  (x-xj) * pow( r , beta - 2.0 ) * ( 1.0 + beta * log(r) ); 
}
//----------------------------------------------------------   
template <typename T> 
T tps_1d_xx(int beta, T x, T xj)
{ 
   T r,tmp1,tmp2;   
   
   r = sqrt( (x-xj) * (x-xj) );   
   
   if( r == 0.0)
       return 0.0;
   else {
      
      tmp1 =  pow(r, beta-4.0);
      tmp2 =  (x-xj) * (x-xj);
      
      return  ( 1.0 + beta*log(r) )  * ( pow(r,  beta-2.0 ) + (beta-2.0) * tmp1 * tmp2) \
               + beta * tmp1 * tmp2; 
   }   
}      
//----------------------------------------------------------   
//              Data  2-D  
//----------------------------------------------------------   
template <typename T> 
T tps_2d_x(int beta, T x, T y, T xj, T yj)
{ 
   T r;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) );   
   
   if( r == 0)
       return 0.0;
   else 
       return  (x-xj) * pow(r,  beta-2.0) * ( 1.0 + beta * log(r) ); 
}
//----------------------------------------------------------   
template <typename T> 
T tps_2d_y(int beta, T x, T y, T xj, T yj)
{ 
   T r;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) );   
   
   if( r == 0)
       return 0.0;
   else 
       return  (y-yj) * pow(r,  beta-2.0) * ( 1.0 + beta * log(r) ); 
}
//----------------------------------------------------------   
template <typename T> 
T tps_2d_xx(int beta, T x, T y, T xj, T yj)
{ 
   T r,tmp1,tmp2;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) );   
   
   if( r == 0.0)
       return 0.0;
   else {
      
      tmp1 =  pow(r, beta-4.0 );
      
      tmp2 =  (x-xj) * (x-xj);
      
      return  ( 1.0 + beta*log(r) )  * ( pow(r, beta-2.0 ) + (beta-2.0) * tmp1 * tmp2) \
      + beta * tmp1 * tmp2; 
   }   
}
//----------------------------------------------------------   
template <typename T> 
T tps_2d_yy(int beta, T x, T y, T xj, T yj)
{ 
   T r,tmp1,tmp2;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) );   
   
   if( r == 0.0)
       return 0.0;
   else {
   	
      tmp1 =  pow(r,  beta-4.0);
      
      tmp2 =  (y-yj) * (y-yj);
      
      return  ( 1.0 + beta*log(r) )  * ( pow(r,  beta-2.0) + (beta-2.0) * tmp1 * tmp2) \
      + beta * tmp1 * tmp2; 
   }   
   
}
//----------------------------------------------------------
template <typename T>
inline T tps_2d_yx(int beta, T x, T y, T xj, T yj)
{
   T  r;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) );   
   
   return pow(r,beta-4.0)*(x-xj)*(y-yj)*( beta + (beta-2.0)*(beta*log(r)+1.0) );
}
//----------------------------------------------------------
template <typename T>
inline T tps_2d_xy(int beta, T x, T y, T xj, T yj)
{
   return tps_2d_yx(beta, x, y, xj, yj);   
}
//----------------------------------------------------------   
//              Data  3-D  
//----------------------------------------------------------   
template <typename T> 
T tps_3d_x(int beta, T x, T y, T xj, T yj, T z, T zj)
{ 
   T r;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj)  + (z-zj) * (z-zj) );   
   
   if( r == 0)
       return 0.0;
   else 
       return  (x-xj) * pow(r,  beta-2.0) * ( 1.0 + beta * log(r) );     
}
//----------------------------------------------------------
template <typename T> 
T tps_3d_y(int beta, T x, T y, T xj, T yj, T z, T zj)
{ 
   T r;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj)  + (z-zj) * (z-zj) );   
   
   if( r == 0)
       return 0.0;
   else 
       return  (y-yj) * pow(r,  beta-2.0) * ( 1.0 + beta * log(r) );     
}
//----------------------------------------------------------
template <typename T> 
T tps_3d_z(int beta, T x, T y, T xj, T yj, T z, T zj)
{ 
   T r;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj)  + (z-zj) * (z-zj) );   
   
   if( r == 0)
       return 0.0;
   else 
       return  (z-zj) * pow(r,  beta-2.0) * ( 1.0 + beta * log(r) );     
}
//----------------------------------------------------------   
template <typename T> 
T tps_2d_xx(int beta, T x, T y, T xj, T yj, T z, T zj)
{ 
   T r,tmp1,tmp2;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj));   
   
   if( r == 0.0)
       return 0.0;
   else {
   	
      tmp1 =  pow(r,  beta-4.0);
      
      tmp2 =  (x-xj) * (x-xj);
      
      return  ( 1.0 + beta*log(r) )  * ( pow(r,  beta-2.0) + (beta-2.0) * tmp1 * tmp2) \
      + beta * tmp1 * tmp2; 
   }   
}
//----------------------------------------------------------   
template <typename T> 
T tps_3d_yy(int beta, T x, T y, T xj, T yj, T z, T zj)
{ 
   T r,tmp1,tmp2;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj));   
   
   if( r == 0.0)
       return 0.0;
   else {
   	
      tmp1 =  pow(r, beta-4.0);
      
      tmp2 =  (y-yj) * (y-yj);
      
      return  ( 1.0 + beta*log(r) )  * ( pow(r, beta-2.0) + (beta-2.0) * tmp1 * tmp2) \
      + beta * tmp1 * tmp2; 
   }   
}
//----------------------------------------------------------   
template <typename T> 
T tps_3d_zz(int beta, T x, T y, T xj, T yj, T z, T zj)
{ 
 T r,tmp1,tmp2;   
   
 r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj));   
   
   if( r == 0.0)
       return 0.0;
   else {
   	
      tmp1 =  pow(r, beta-4.0);
      
      tmp2 =  (z-zj) * (z-zj);
      
      return  ( 1.0 + beta*log(r) )  * ( pow(r,  beta-2.0) + (beta-2.0) * tmp1 * tmp2) \
      + beta * tmp1 * tmp2; 
   }   
}

} // RBF namespace

#endif //_RBF_TPS_DERIVATIVES_HPP_




