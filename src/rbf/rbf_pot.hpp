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

#ifndef _RBF_POT_CPP_
#define _RBF_POT_CPP_
#include "rbf_pot.h"

namespace rbf{

//----------------------------------------------------------
template <typename T> 
POT<T>::POT(void)
{
   beta            = 1; 
   degree          = get_min_degree_pol();
   initialized     = false;
   initialized_pol = false;
}
//----------------------------------------------------------------
template <class T>
std::ostream& operator<<(std::ostream &s,  POT<T> &rbf)
{
    s <<"RBF info: "<<"\n";
    s <<" name           : " << rbf.name()<< "\n";  
    s <<" beta factor    : " << rbf.get_beta()<< "\n";
    s <<" degree pol     : " << rbf.get_degree_pol()<< "\n";
    s <<" min degree pol : " << rbf.get_min_degree_pol()<< "\n"; 
    return s;
}
//----------------------------------------------------------
template <typename T> 
void POT<T>::set_beta(int beta_factor)
{
   //beta = 1, 3, 5, 7, 9, 11, ...
   
   if(beta_factor<1)
   {
      fprintf(stderr,"\nERROR in POT-set_beta, beta = %d, must be a positive integer  1,3,5,... \n\n",beta_factor);
      exit(1);
   }   
   
   if(!rbf_is_impar(beta_factor))
   {
      fprintf(stderr,"\nERROR: in POT-set_beta, beta = %d, must be  1,3,5,7,... \n\n",beta_factor);
      exit(1);
   }      
   
   beta        = beta_factor; 
   
   if(initialized_pol)
         degree      = degree ;
    else
         degree      =  get_min_degree_pol();
        
   initialized = true;
}
//----------------------------------------------------------
template <typename T>
int POT<T>::get_beta(void)
{
  return int(beta);   
}
//----------------------------------------------------------
template<typename T>
int POT<T>::get_degree_pol(void)
{
  return degree;   
}
//----------------------------------------------------------
template <typename T>
void  POT<T>::set_degree_pol(int degree_factor)
{

  if(degree_factor<0){
    fprintf(stderr,"\n!! ERROR in POT-set_degree_pol: the degree = %d must be >= 0.\n\n",degree_factor);
    fflush(stderr);
    exit(1);
  }
  
   //First, validate that degre_factor >= min_degree_pol requiered
   if( degree_factor < get_min_degree_pol() )
   {
      fprintf(stderr,"\nWARNING: in POT kernel the degree = %d must be at least %d.\n\n",degree_factor, get_min_degree_pol());
   }   
   
   degree = degree_factor;   
   
   initialized_pol = true;
}
//----------------------------------------------------------
template<typename T>
int POT<T>::get_min_degree_pol(void)
{ 
  return int( ceil(  T(beta) / 2.0) );   
}
//----------------------------------------------------------
template <typename T> 
inline T POT<T>::eval(T x,T xj, T c)
{
   T r;
   
   r = sqrt( (x-xj)*(x-xj) );      
   
 if(initialized==false)
   return   r; 
 else
   return rbf_pow( r , beta ); 
  
}
//----------------------------------------------------------
template <typename T> 
inline T POT<T>::eval(T x, T y, T xj, T yj, T c)
{
   T r;
   
   r = sqrt( (x-xj)*(x-xj) + (y-yj)*(y-yj) );      
   
 if(initialized==false)
   return   r ; 
 else
   return rbf_pow( r  , beta ); 

}
//----------------------------------------------------------
template <typename T> 
inline T POT<T>::eval(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r;
   
   r = sqrt( (x-xj)*(x-xj) + (y-yj)*(y-yj) +(z-zj)*(z-zj) );   
   
   if(initialized==false)
      return   r ; 
   else 
      return rbf_pow(  r , beta ); 
}
//----------------------------------------------------------
template <typename T> 
T POT<T>::eval(const T *x, const T *xj, int dim, T c)
{ 
   T r;
   
   r = sqrt(rbf_norm_square(x,xj,dim));
   
   if(initialized==false)
      return   r ; 
   else
      return rbf_pow(  r ,  beta ); 
}
//----------------------------------------------------------   
//              Data  1-D  
//----------------------------------------------------------
template <typename T>
inline T POT<T>::dx(T x, T xj, T c)
{  
  T r = sqrt( (x-xj)*(x-xj) );
  
   return (x-xj) * beta * rbf_pow( r , beta-2 ) ;
}
//----------------------------------------------------------
template <typename T>
inline T POT<T>::dxx(T x, T xj, T c)
{  
  T r = sqrt( (x-xj)*(x-xj) );
  
  
   return  beta * rbf_pow( r , beta-2 ) + beta * (x-xj) * (x-xj) * (beta-2.0) * rbf_pow( r , beta-4.0);
}
//----------------------------------------------------------   
//              Data  2-D  
//----------------------------------------------------------
 template <typename T>
inline T POT<T>::dx(T x, T y, T xj, T yj, T c)
{
   T   r;
   
    r = sqrt( ( x-xj) * (x-xj) + (y-yj) * (y-yj) );
  
   return (x-xj) * beta * rbf_pow( r , beta-2 ) ;

}
//----------------------------------------------------------
template <typename T>
inline T POT<T>::dy(T x, T y, T xj, T yj, T c)
{
   T   r;
   
    r = sqrt( ( x-xj)*(x-xj) + (y-yj)*(y-yj) );
  
   return (y-yj) * beta * rbf_pow( r , beta-2 ) ;
}
//----------------------------------------------------------
template <typename T>
inline T POT<T>::dxx(T x, T y, T xj, T yj, T c)
{
   T   r;
   
    r = sqrt( ( x-xj)*(x-xj) + (y-yj)*(y-yj) );

   return  beta * rbf_pow( r , beta-2 ) + beta * (x-xj) * (x-xj) * (beta-2) * rbf_pow( r , beta-4);  
}
//----------------------------------------------------------
template <typename T>
inline T POT<T>::dyy(T x, T y, T xj, T yj, T c)
{
   T   r;
   
    r = sqrt( ( x-xj)*(x-xj) + (y-yj)*(y-yj) );

   return  beta * rbf_pow( r , beta-2 ) + beta * (y-yj) * (y-yj) * (beta-2) * rbf_pow( r , beta-4);  
}
//----------------------------------------------------------
template <typename T>
inline T POT<T>::dyx(T x, T y, T xj, T yj, T c)
{
   T  r;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) );   

   return beta * (beta-1.0) * (x-xj) * (y-yj) * pow(r,beta-3.0) ;
}	
//----------------------------------------------------------
template <typename T>
inline T POT<T>::dxy(T x, T y, T xj, T yj, T c)
{
   T  r;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) );   

   return beta * (beta-1.0) * (x-xj) * (y-yj) * pow(r,beta-3.0) ;
}
//----------------------------------------------------------   
//              Data  3-D  
//----------------------------------------------------------   
template <typename T> 
T POT<T>::dx(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T   r;
   
   r = sqrt( ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj) );
  
   return (x-xj) * beta * pow( r , beta-2.0 ) ;
}
//----------------------------------------------------------
template <typename T> 
T POT<T>::dy(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T   r;
   
   r = sqrt( ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj) );
  
   return (y-yj) * beta * pow( r , beta-2.0 ) ;
}
//----------------------------------------------------------
template <typename T> 
T POT<T>::dz(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T   r;
   
   r = sqrt( ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj) );
  
   return (z-zj) * beta * pow( r , beta-2.0 ) ;
}
//----------------------------------------------------------
template <typename T> 
T POT<T>::dxx(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T   r, bm2;
   
   r   = sqrt( ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj) );
   
   bm2 = beta - 2.0;
   
   return  beta * pow( r , bm2 ) + beta * bm2 * pow( r , beta-4.0 ) * (x-xj) * (x-xj) ;  
}
//----------------------------------------------------------
template <typename T> 
T POT<T>::dyy(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T   r, bm2;
   
   r   = sqrt( ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj) );
   
   bm2 = beta - 2.0;
   
   return  beta * pow( r , bm2 ) + beta * bm2 * pow( r , beta-4.0 ) * (y-yj) * (y-yj) ;  
}
//----------------------------------------------------------
template <typename T> 
T POT<T>::dzz(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T   r, bm2;
   
   r   = sqrt( ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj) );
   
   bm2 = beta - 2.0;
   
   return  beta * pow( r , bm2 ) + beta * bm2 * pow( r , beta-4.0 ) * (z-zj) * (z-zj) ;  
}



} // RBF namespace

#endif // _RBF_POT_CPP_

