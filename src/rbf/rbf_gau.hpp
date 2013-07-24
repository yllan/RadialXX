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


#ifndef _RBF_GAU_CPP_
#define _RBF_GAU_CPP_

#include "rbf_gau.h"

namespace rbf{

//----------------------------------------------------------
template <typename T> 
GAU<T>::GAU(void)
{
   beta            = 0; 
   degree          = get_min_degree_pol();
   initialized     = false;
   initialized_pol = false;
}
//----------------------------------------------------------------
template <class T>
std::ostream& operator<<(std::ostream &s,  GAU<T> &rbf)
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
void GAU<T>::set_beta(int beta_factor)
{
   fprintf(stderr,"WARNING: in GAU kernel, does not requiere the beta factor.\n\n");
   fflush(stderr);
   
   beta        = beta_factor; 
   
   if(initialized_pol)
    degree      = degree ;
    else
    degree      = get_min_degree_pol();
         
  initialized = true;
}
//----------------------------------------------------------
template <typename T>
int GAU<T>::get_beta(void)
{
   return int(beta);   
}
//----------------------------------------------------------
template<typename T>
int GAU<T>::get_degree_pol(void)
{
   return degree;   
}
//----------------------------------------------------------
template <typename T>
void  GAU<T>::set_degree_pol(int degree_factor)
{

  if(degree_factor<0){
    fprintf(stderr,"\n!! ERROR in GAU-set_degree_pol: the degree = %d must be >= 0.\n\n",degree_factor);
    fflush(stderr);
    exit(1);
  }

   //First, validate that degre_factor >= min_degree_pol requiered
   if( degree_factor < get_min_degree_pol() )
   {
      fprintf(stderr,"WARNING: in GAU kernel the degree = %d must be at least %d.\n\n",
                     degree_factor,   get_min_degree_pol());
   }   
   
   degree = degree_factor;   
   
   initialized_pol = true;
}
//----------------------------------------------------------
template<typename T>
int GAU<T>::get_min_degree_pol(void)
{
  return 0;   
}
//----------------------------------------------------------
template <typename T> 
T GAU<T>::eval(T x,T xj, T c)
{
   T r2;
   
   r2 = (x-xj)*(x-xj);      
   
    return exp( -r2*c*c); 
  
}
//----------------------------------------------------------
template <typename T> 
T GAU<T>::eval(T x, T y, T xj, T yj, T c)
{
   T r2;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj);      
   
    return exp(-r2*c*c);
}
//----------------------------------------------------------
template <typename T> 
T GAU<T>::eval(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj) +(z-zj)*(z-zj);   

   return exp(-r2*c*c);
}
//----------------------------------------------------------
template <typename T> 
T GAU<T>::eval(const T *x,const T *xj, int dim, T c)
{ 
   T r;
      
   r = rbf_norm_square(x,xj,dim);
   
    return exp(-r*c*c); 
}
//----------------------------------------------------------   
//              Data  1-D  
//----------------------------------------------------------   
template <typename T>
inline T GAU<T>::dx(T x, T xj, T c)
{  
   T r2;
      
   r2 =   (x - xj) * (x - xj);
   
   return -2.0 * c*c * (x - xj) * exp( -r2 * c*c );
}
//----------------------------------------------------------
template <typename T>
inline T GAU<T>::dxx(T x, T xj, T c)
{  
   T r2;
      
   r2 =   (x - xj) * (x - xj);
   
   return 2.0 * c*c * exp( -r2 * c*c ) * ( 2.0 * c*c * (x - xj) * (x - xj) - 1.0 );
}
//----------------------------------------------------------   
//              Data  2-D  
//----------------------------------------------------------
 template <typename T>
inline T GAU<T>::dx(T x, T y, T xj, T yj, T c)
{
   T r2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj);
   
   return -2.0 * c*c * (x - xj) * exp( -r2 * c*c );
}
//----------------------------------------------------------
template <typename T>
inline T GAU<T>::dy(T x, T y, T xj, T yj, T c)
{
   T r2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj);
   
   return -2.0 * c*c * (y - yj) * exp( -r2 * c*c );
}
//----------------------------------------------------------
template <typename T>
inline T GAU<T>::dxx(T x, T y, T xj, T yj, T c)
{
   T r2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj);
   
   
   return 2.0 * c*c * exp( -r2 * c*c ) * ( 2.0 * c*c * (x - xj) * (x - xj) - 1.0 );
   
}
//----------------------------------------------------------
template <typename T>
inline T GAU<T>::dyy(T x, T y, T xj, T yj, T c)
{
   T r2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj);
   
   return 2.0 * c*c * exp( -r2 * c*c ) * ( 2.0 * c*c * (y - yj) * (y - yj) - 1.0 ); 
}
//----------------------------------------------------------
template <typename T>
inline T GAU<T>::dyx(T x, T y, T xj, T yj, T c)
{
   T  r2, c2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj);
   
   c2 = c * c;

   return (-2.0 * c2 * (y-yj) )*(-2.0 * c2 * (x-xj) )*exp(-r2 * c2);
}
//----------------------------------------------------------
template <typename T>
inline T GAU<T>::dxy(T x, T y, T xj, T yj, T c)
{
   T  r2, c2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj);
   
   c2 = c * c;

   return (-2.0 * c2 * (y-yj) )*(-2.0 * c2 * (x-xj) )*exp(-r2 * c2);
}
//----------------------------------------------------------   
//              Data  3-D  
//----------------------------------------------------------   
template <typename T> 
T GAU<T>::dx(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj);
   
   return -2.0 * c*c * (x - xj) * exp( -r2 * c*c );  
}
//----------------------------------------------------------
template <typename T> 
T GAU<T>::dy(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj);
   
   return -2.0 * c*c * (y - yj) * exp( -r2 * c*c );  
}
//----------------------------------------------------------
template <typename T> 
T GAU<T>::dz(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj);
   
   return -2.0 * c*c * (z - zj) * exp( -r2 * c*c );  
}
//----------------------------------------------------------
template <typename T> 
T GAU<T>::dxx(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2,c2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj);
   
   c2 = c * c;
   
   return 2.0 * c2 * exp( -r2 * c2 ) * ( 2.0 * c2 * (x - xj) * (x - xj) - 1.0 );
}
//----------------------------------------------------------
template <typename T> 
T GAU<T>::dyy(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2,c2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj);
   
   c2 = c * c;
   
   return 2.0 * c2 * exp( -r2 * c2 ) * ( 2.0 * c2 * (y - yj) * (y - yj) - 1.0 );
}
//----------------------------------------------------------
template <typename T> 
T GAU<T>::dzz(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2,c2;
      
   r2 =   (x - xj) * (x - xj) + (y - yj) * (y - yj) + (z - zj) * (z - zj);
   
   c2 = c * c;
   
   return 2.0 * c2 * exp( -r2 * c2 ) * ( 2.0 * c2 * (z - zj) * (z - zj) - 1.0 );
}


} // RBF namespace

#endif // _RBF_GAU_CPP_



