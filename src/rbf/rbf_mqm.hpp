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

#ifndef _RBF_MQM_HPP_
#define _RBF_MQM_HPP_

#include "rbf_mqm.h"

namespace rbf{

//----------------------------------------------------------
template <typename T> 
MQM<T>::MQM(void)
{
  beta            = 1; 
  degree          = get_min_degree_pol();
  initialized     = false;
  initialized_pol = false;
}
//----------------------------------------------------------------
template <class T>
std::ostream& operator<<(std::ostream &s,  MQM<T> &rbf)
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
void MQM<T>::set_beta(int beta_factor)
{
   // beta = 1, 3, 5, 7, 9, 11, ...
   // phi(r) =  (r^2 + c^2)^(beta/2)   
   
  if(beta_factor<1)
  {
   fprintf(stderr,"ERROR: in MQ kernel, beta = %d debe ser un numero positivo. \n\n",beta_factor);
   exit(1);
  }   
   
  if(!rbf_is_impar(beta_factor))
  {
   fprintf(stderr,"ERROR: in MQ kernel, beta = %d debe ser un numero impar. \n\n",beta_factor);
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
int MQM<T>::get_beta(void)
{
  return int(beta);   
}
//----------------------------------------------------------
template<typename T>
int MQM<T>::get_degree_pol(void)
{
  return degree;   
}

//----------------------------------------------------------
template <typename T>
void  MQM<T>::set_degree_pol(int degree_factor)
{

  if(degree_factor<0){
    fprintf(stderr,"\n!! ERROR in MQ-set_degree_pol: the degree = %d must be >= 0.\n\n",degree_factor);
    fflush(stderr);
    exit(1);
  }
  
  //First, validate that degre_factor >= min_degree_pol requiered
  if( degree_factor < get_min_degree_pol() )
  {
    fprintf(stderr,"\nWARNING: in MQ kernel the degree assigned = %d must be at least %d.\n\n",degree_factor, get_min_degree_pol());
  }   
   
  degree = degree_factor;   
  
  initialized_pol = true;
    
}
//----------------------------------------------------------
template<typename T>
int MQM<T>::get_min_degree_pol(void)
{
  // incorrect  return int(beta);   //WARNING it is not vaidated with the theory
  return int( ceil(  beta / 2.0) );   
}
//----------------------------------------------------------
template <typename T> 
inline T MQM<T>::eval(T x,T xj, T c)
{
   T c2=c*c;   
   T r2;
   
   r2 = (x-xj)*(x-xj);      
   
   if(beta==1)
      return sqrt( 1.0  + c2 * r2  ); 
   else
      return pow( 1.0 + c2 * r2 ,   beta / 2.0 ); 
  
}
//----------------------------------------------------------
template <typename T> 
inline T MQM<T>::eval(T x, T y, T xj, T yj, T c)
{
   T c2=c*c;
   T r2;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj);      
   
   if(beta==1)
      return sqrt( 1.0 + c2 * r2 ); 
   else
      return pow( 1.0 + c2 * r2 ,  beta / 2.0 ); 

}
//----------------------------------------------------------
template <typename T> 
inline T MQM<T>::eval(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T c2=c*c;
   T r2;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj) +(z-zj)*(z-zj);
   
   if(beta==1)
     return sqrt( 1.0 + c2 * r2 ); 
   else
     return pow(  1.0 + c2 * r2 ,   beta / 2.0 ); 
}
//----------------------------------------------------------
template <typename T> 
inline T MQM<T>::eval(const T *x,const T *xj, int dim, T c)
{ 
   T c2=c*c;   
   T r2;

   r2 = rbf_norm_square(x,xj,dim);

   if(beta==1)
       return sqrt( 1.0 + c2 * r2 ); 
   else
       return pow( 1.0 + c2 * r2 ,  beta / 2.0); 
}



} // RBF namespace

#endif // _RBF_MQM_HPP_
