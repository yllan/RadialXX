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


#ifndef _RBF_TPS_CPP_
#define _RBF_TPS_CPP_

#include "rbf_tps.h"

namespace rbf{
   
//----------------------------------------------------------
template <typename T> 
TPS<T>::TPS(void)
{
   beta            = 2; 
   degree          = get_min_degree_pol();
   initialized     = false;
   initialized_pol = false;
}
//----------------------------------------------------------------
template <class T>
std::ostream& operator<<(std::ostream &s,  TPS<T> &rbf)
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
void TPS<T>::set_beta(int beta_factor)
{
   //beta = 2, 4, 6, 8, 10, ...
   if(beta_factor<2)
   {
       fprintf(stderr,"\nERROR: in TPS kernel, beta = %d debe ser un numero positivo par mayor o igual a 2. \n\n",beta_factor);
       fflush(stderr);
       exit(1);
   }      
   
   if(rbf_is_impar(beta_factor))
   {
         fprintf(stderr,"\nERROR: beta = %d debe ser un numero par. \n\n",beta_factor);
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
int TPS<T>::get_beta(void)
{
   return int(beta);   
}
//----------------------------------------------------------
template<typename T>
int TPS<T>::get_degree_pol(void)
{
   return degree;   
}
//----------------------------------------------------------
template <typename T>
void  TPS<T>::set_degree_pol(int degree_factor)
{
   if(degree_factor<0){
      fprintf(stderr,"\n!! ERROR in TPS-set_degree_pol: the degree = %d must be >= 0.\n\n",degree_factor);
      fflush(stderr);
      exit(1);
    }


//First, validate that degre_factor >= min_degree_pol requiered
   if( degree_factor < get_min_degree_pol() )
   {
      fprintf(stderr,"\nWARNING: in TPS kernel the degree = %d must be at least %d.\n\n",degree_factor, get_min_degree_pol());
      fflush(stderr);
   }   
   
   degree = degree_factor;   
   
   initialized_pol = true;   
}
//----------------------------------------------------------
template<typename T>
int TPS<T>::get_min_degree_pol( void)
{
   return int( 1.0 + beta/2.0 );   
}

//----------------------------------------------------------
template <typename T> 
T TPS<T>::eval(T x,T xj, T c)
{
   T r;
   
   r = sqrt( (x-xj)*(x-xj) );      
   
   if(r==0)
      return 0;
   else      
      return rbf_pow( r , int(beta) )*log(r); 
}
//----------------------------------------------------------
template <typename T> 
T TPS<T>::eval(T x, T y, T xj, T yj, T c)
{
   T r;

   
   r = sqrt( (x-xj)*(x-xj) + (y-yj)*(y-yj) );      
   
   if(r==0)
      return 0;
   else
      return rbf_pow( r  , int(beta) )*log(r); 
}
//----------------------------------------------------------
template <typename T> 
T TPS<T>::eval(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r;
   
   r = sqrt( (x-xj)*(x-xj) + (y-yj)*(y-yj) +(z-zj)*(z-zj) );   
   
   if(r==0)
      return 0;
   else
      return rbf_pow(  r , int(beta) )*log(r); 
}
//----------------------------------------------------------
template <typename T> 
T TPS<T>::eval(const T *x, const T *xj, int dim, T c)
{ 
   T r;
   
   r = sqrt( rbf_norm_square(x,xj,dim));   
   
   if(r==0)
      return 0;
     else 
      return rbf_pow(  r , int(beta) )*log(r); 
}
//----------------------------------------------------------   
//              Data  1-D  
//----------------------------------------------------------   
template <typename T> 
T TPS<T>::dx(T x, T xj,  T c)
{ 
   T r;   
   
   r = sqrt( (x-xj) * (x-xj) );   
   
   if( r == 0)
       return 0.0;
   else 
       return  (x-xj) * pow(r,  beta - 2.0) * ( 1.0 + beta * log(r) ); 
}
//----------------------------------------------------------   
template <typename T> 
T TPS<T>::dxx(T x, T xj,  T c)
{ 
   T r,tmp1,tmp2;   
   
   r = sqrt( (x-xj) * (x-xj) );   
   
   if( r == 0.0)
       return 0.0;
   else {
      tmp1 =  pow(r, (T)(beta-4.0));
      tmp2 =  (x-xj) * (x-xj);
       return  ( 1.0 + beta*log(r) )  * ( pow(r, (T)(beta-2.0)) + (beta-2.0) * tmp1 * tmp2) \
               + beta * tmp1 * tmp2; 
   }   
}      
//----------------------------------------------------------   
//              Data  2-D  
//----------------------------------------------------------   
template <typename T> 
T TPS<T>::dx(T x, T y, T xj, T yj, T c)
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
T TPS<T>::dy(T x, T y, T xj, T yj, T c)
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
T TPS<T>::dxx(T x, T y, T xj, T yj, T c)
{ 
   T r,tmp1,tmp2;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) );   
   
   if( r == 0.0)
       return 0.0;
   else {
      tmp1 =  pow(r, (T)(beta-4.0));
      tmp2 =  (x-xj) * (x-xj);
       return  ( 1.0 + beta*log(r) )  * ( pow(r, (T)(beta-2.0)) + (beta-2.0) * tmp1 * tmp2) \
      + beta * tmp1 * tmp2; 
   }   
}
//----------------------------------------------------------   
template <typename T> 
T TPS<T>::dyy(T x, T y, T xj, T yj, T c)
{ 
   T r,tmp1,tmp2;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) );   
   
   if( r == 0.0)
       return 0.0;
   else {
      tmp1 =  pow(r, (T)(beta-4.0));
      tmp2 =  (y-yj) * (y-yj);
       return  ( 1.0 + beta*log(r) )  * ( pow(r,  beta-2.0) + (beta-2.0) * tmp1 * tmp2) \
      + beta * tmp1 * tmp2; 
   }   
   
}
//----------------------------------------------------------   
template <typename T> 
T TPS<T>::dyx(T x, T y, T xj, T yj, T c)
{ 
   T  r;   
   
   r = sqrt( (x-xj) * (x-xj) + (y-yj) * (y-yj) );   
   
   return pow(r,beta-4.0)*(x-xj)*(y-yj)*( beta + (beta-2.0)*(beta*log(r)+1.0) );
}	
//----------------------------------------------------------   
template <typename T> 
T TPS<T>::dxy(T x, T y, T xj, T yj, T c)
{ 
   return this->dyx(x, y, xj, yj, c);   
}
//----------------------------------------------------------   
//              Data  3-D  
//----------------------------------------------------------   
template <typename T> 
T TPS<T>::dx(T x, T y, T z, T xj, T yj, T zj, T c)
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
T TPS<T>::dy(T x, T y, T z, T xj, T yj, T zj, T c)
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
T TPS<T>::dz(T x, T y, T z, T xj, T yj, T zj, T c)
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
T TPS<T>::dxx(T x, T y, T z, T xj, T yj, T zj, T c)
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
T TPS<T>::dyy(T x, T y, T z, T xj, T yj, T zj, T c)
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
T TPS<T>::dzz(T x, T y, T z, T xj, T yj, T zj, T c)
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

#endif //_RBF_TPS_CPP_


