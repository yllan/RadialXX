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

// Perform the evaluation 1d, 2d, 3d and nd of the generalized radial basis functions
//
//  multiquadric
//  inverse multiquadric
//  gaussian
//  thin-plate splines
//  power splines
//

#ifndef _RBF_EVAL_HPP_
#define _RBF_EVAL_HPP_


namespace rbf{
	
//----------------------------------------------------------   
//             Thin-Plate Splines  
//---------------------------------------------------------- 
template <typename T> 
T tps_1d(int beta , T x,T xj)
{
   T r;
   
   r = sqrt( (x-xj)*(x-xj) );      
   
   if(r==0)
      return 0;
   else      
      return rbf_pow( r , beta )*log(r); 
}
//
//----------------------------------------------------------
//
template <typename T> 
T tps_2d(int beta, T x, T y, T xj, T yj, T c)
{
   T r;
   
   r = sqrt( (x-xj)*(x-xj) + (y-yj)*(y-yj) );      
   
   if(r==0)
      return 0;
   else
      return rbf_pow( r  , beta )*log(r); 
}
//
//----------------------------------------------------------
//
template <typename T> 
T tps_3d(int beta, T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r;
   
   r = sqrt( (x-xj)*(x-xj) + (y-yj)*(y-yj) + (z-zj)*(z-zj) );   
   
   if(r==0)
      return 0;
   else
      return rbf_pow(  r , beta )*log(r); 
}
//
//----------------------------------------------------------
//
template <typename T> 
T tps_nd(int beta, const T *x, const T *xj, int dim, T c)
{ 
   T r;
   
   r = sqrt( rbf_norm_square(x,xj,dim));   
   
   if(r==0)
      return 0;
    else 
      return rbf_pow(  r , beta )*log(r); 
}
//----------------------------------------------------------   
//             Power Splines  
//---------------------------------------------------------- 
template <typename T> 
T pot_1d(int beta , T x,T xj)
{
   T r;
   
   r = sqrt( (x-xj)*(x-xj) );      
       
   return rbf_pow( r ,  beta ); 
}
//
//----------------------------------------------------------
//
template <typename T> 
T pot_2d(int beta, T x, T y, T xj, T yj, T c)
{
   T r;
   
   r = sqrt( (x-xj)*(x-xj) + (y-yj)*(y-yj) );      
   
   return rbf_pow( r ,  beta );  
}
//
//----------------------------------------------------------
//
template <typename T> 
T pot_3d(int beta, T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r;
   
   r = sqrt( (x-xj)*(x-xj) + (y-yj)*(y-yj) +(z-zj)*(z-zj) );   
   
   return rbf_pow( r ,  beta );  
}
//
//----------------------------------------------------------
//
template <typename T> 
T pot_nd(int beta, const T *x, const T *xj, int dim, T c)
{ 
   T r;
   
   r = sqrt( rbf_norm_square(x,xj,dim) );   
   
   return rbf_pow( r ,  beta ); 
}
//----------------------------------------------------------   
//              Gaussian  Kernel  
//----------------------------------------------------------   
template <typename T>
inline T gau_1d(T x, T xj, T c)
{  
   T r2;
      
   r2 = (x - xj) * (x - xj);
   
   return  exp( -r2 * c*c );
}
//
//----------------------------------------------------------
//
template <typename T>
inline T gau_2d(T x, T y, T xj, T yj, T c)
{  
   T r2;
      
   r2 = (x - xj) * (x - xj) + (y - yj) * (y - yj);    
   
   return  exp( -r2 * c*c );
}
//
//----------------------------------------------------------
//
template <typename T>
inline T gau_3d(T x, T y, T z, T xj, T yj, T zj, T c)
{  
   T r2;
      
   r2 =  (x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj);
   
   return  exp( -r2 * c*c );
}
//
//----------------------------------------------------------
//
template <typename T>
inline T gau_nd(const T *x,const T *xj, int dim, T c)
{  
   T r;
      
   r = rbf_norm_square(x,xj,dim);
   
   return exp(-r*c*c); 
}
//----------------------------------------------------------   
//             Multiquadric 
//---------------------------------------------------------- 
template <typename T> 
inline T mq_1d(int beta, T x,T xj, T c)
{
   T r2;
   
   r2 = (x-xj)*(x-xj);      
   
   if(beta==1)
      return sqrt( r2 + c*c ); 
   else
      return pow( r2 + c*c,   beta / 2.0 ); 
  
}
//----------------------------------------------------------
template <typename T> 
inline T mq_2d(int beta, T x, T y, T xj, T yj, T c)
{
   T r2;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj);      
   
   if(beta==1)
      return sqrt( r2  + c*c); 
   else
      return pow( r2  + c*c ,  beta / 2.0 ); 

}
//----------------------------------------------------------
template <typename T> 
inline T mq_3d(int beta, T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj) +(z-zj)*(z-zj);
   
   if(beta==1)
     return sqrt( r2 + c*c); 
   else
     return pow(  r2  + c*c ,   beta / 2.0 ); 
}
//
//----------------------------------------------------------
//
template <typename T> 
inline T mq_nd(int beta, const T *x,const T *xj, int dim, T c)
{ 
   T r2;

   r2 = rbf_norm_square(x,xj,dim);

   if(beta==1)
      return sqrt( r2 + c*c); 
   else
      return pow(  r2  + c*c ,  beta / 2.0); 
}
//----------------------------------------------------------   
//            Inverse Multiquadric 
//---------------------------------------------------------- 
template <typename T> 
inline T imq_1d(int beta, T x,T xj, T c)
{
   T r2;
   
   r2 = (x-xj)*(x-xj);      
   
   if(beta==1)
      return 1.0/sqrt( r2 + c*c ); 
   else
      return pow( r2 + c*c, -beta/2.0); 
  
}
//----------------------------------------------------------
template <typename T> 
inline T imq_2d(int beta, T x, T y, T xj, T yj, T c)
{
   T r2;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj);      
   
   if(beta==1)
      return 1.0/sqrt( r2  + c*c); 
   else
      return pow( r2  + c*c , -beta/2.0); 

}
//----------------------------------------------------------
template <typename T> 
inline T imq_3d(int beta, T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj) +(z-zj)*(z-zj);   
   
   if(beta==1)
      return 1.0/sqrt( r2 + c*c); 
   else
      return pow(  r2  + c*c , -beta/2.0 ); 
}
//
//----------------------------------------------------------
//
template <typename T> 
inline T imq_nd(int beta, const T *x,const T *xj, int dim, T c)
{ 
   T r2;
   
   r2 = rbf_norm_square(x,xj,dim);   
   
   if(beta==1)
      return 1.0/sqrt( r2 + c*c); 
   else
      return pow(  r2 + c*c , -beta/2.0); 
}


} // RBF namespace

#endif // _RBF_EVAL_HPP_





