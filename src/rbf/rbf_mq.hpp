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

#ifndef _RBF_MQ_CPP_
#define _RBF_MQ_CPP_

#include "rbf_mq.h"

namespace rbf{

//----------------------------------------------------------
template <typename T> 
MQ<T>::MQ(void)
{
  beta            = 1; 
  degree          = get_min_degree_pol();
  initialized     = false;
  initialized_pol = false;
}
//----------------------------------------------------------------
template <class T>
std::ostream& operator<<(std::ostream &s,  MQ<T> &rbf)
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
void MQ<T>::set_beta(int beta_factor)
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
int MQ<T>::get_beta(void)
{
  return int(beta);   
}
//----------------------------------------------------------
template<typename T>
int MQ<T>::get_degree_pol(void)
{
  return degree;   
}

//----------------------------------------------------------
template <typename T>
void  MQ<T>::set_degree_pol(int degree_factor)
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
int MQ<T>::get_min_degree_pol(void)
{
  // incorrect  return int(beta);   //WARNING it is not vaidated with the theory
  return int( ceil(  beta / 2.0) );   
}
//----------------------------------------------------------
template <typename T> 
inline T MQ<T>::eval(T x,T xj, T c)
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
inline T MQ<T>::eval(T x, T y, T xj, T yj, T c)
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
inline T MQ<T>::eval(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
  T r2;
   
  r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj) +(z-zj)*(z-zj);
   
  if(beta==1)
   return sqrt( r2 + c*c); 
  else
   return pow(  r2  + c*c ,   beta / 2.0 ); 
}
//----------------------------------------------------------
template <typename T> 
inline T MQ<T>::eval(const T *x,const T *xj, int dim, T c)
{ 
  T r2;

  r2 = rbf_norm_square(x,xj,dim);

  if(beta==1)
      return sqrt( r2 + c*c); 
  else
      return pow(  r2  + c*c ,  beta / 2.0); 
}
//----------------------------------------------------------   
//              Data  1-D  
//----------------------------------------------------------
template <typename T>
inline T MQ<T>::dx(T x, T xj, T c)
{  
   T r2;
   
   r2 = (x-xj) * (x-xj);
   
   if(beta==1)
      return (x-xj) / sqrt( r2 + c*c );
   else
      return  beta * (x-xj) * pow( r2 + c*c,  beta / 2.0 - 1.0  );   
}
//----------------------------------------------------------
template <typename T>
inline T MQ<T>::dxx(T x, T xj, T c)
{  
   T r2,c2;
   
   c2 = c*c;
   r2 = (x-xj)*(x-xj);
   
   if(beta==1)
      return  (c*c) / ( sqrt(r2 + c*c) * (r2 + c*c) ); 
   else   
       return  beta * pow( r2 + c2,  beta/2.0 - 1.0) * ( 1.0 + ( ( beta-2.0) * r2 )/(r2 + c2) );   
}
//----------------------------------------------------------   
//              Data  2-D  
//----------------------------------------------------------
 template <typename T>
inline T MQ<T>::dx(T x, T y, T xj, T yj, T c)
{
   T r2;
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj); 
   
   if(beta==1)
      return (x-xj) / sqrt( r2 + c*c );
   else   
   
   return   beta * (x-xj) * pow( r2 + c*c,   beta/2.0 - 1.0);   
   
}
//----------------------------------------------------------
template <typename T>
inline T MQ<T>::dy(T x, T y, T xj, T yj, T c)
{
   T r2;
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj); 
   
   
   if(beta==1)
       return (y-yj) / sqrt( r2 + c*c );
   else
       return  beta * (y-yj) * pow( r2 + c*c,  beta / 2.0 - 1.0);      
   
}
//----------------------------------------------------------
template <typename T>
inline T MQ<T>::dxx(T x, T y, T xj, T yj, T c)
{
   T r2,c2;
   T den1,den2;
   c2 = c*c;
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj);
   
   
   if(beta==1){
      //return (1.0/(sqrt(r2 + c2) )) * (1.0  - ((x-xj)*(x-xj))/(r2+c2) ); 
        
        den1=std::pow( r2+c*c,(T)(3.0/2.0));
        den2=sqrt(r2+c*c);
        return  -((x-xj)*(x-xj))/den1 + 1.0/den2;
        
   }    
  else
       return   beta *  pow( r2 + c2,  beta/2.0 - 1.0 ) * ( 1.0 + ( ( beta-2.0) * ( (x-xj)*(x-xj) ) )/(r2 + c2) );
   
}
//----------------------------------------------------------
template <typename T>
inline T MQ<T>::dyy(T x, T y, T xj, T yj, T c)
{
   T r2,c2;
   T den1,den2;
   
   c2 = c*c;
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj);
   
   
   if(beta==1){
      //return (1.0/(sqrt(r2 + c2) )) * (1.0  - ((y-xj)*(y-xj))/(r2+c2) ); 
        
        den1=std::pow(r2+c*c,(T)(3.0/2.0));
        den2=sqrt(r2+c*c);
        return  -((y-yj)*(y-yj))/den1 + 1.0/den2;
        
   }   
   else   
      return   beta * pow( r2 + c2,  beta/2.0 - 1.0 ) * ( 1.0 + ( (  beta - 2.0) * ( (x-xj)*(x-xj) ) )/(r2 + c2) );   
   
}
//----------------------------------------------------------
template <typename T>
inline T MQ<T>::dyx(T x, T y, T xj, T yj, T c)
{
   T  r2;   
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj) ;   
   
   return  beta * (beta-2.0) * (x-xj) * (y-yj) * pow( r2 + c*c,  beta / 2.0 - 2.0);    
}	
//----------------------------------------------------------
template <typename T>
inline T MQ<T>::dxy(T x, T y, T xj, T yj, T c)
{
	return this->dyx(x,y,xj,yj,c);
}	
//----------------------------------------------------------   
//              Data  3-D  
//----------------------------------------------------------   
template <typename T> 
T MQ<T>::dx(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2;
   
   r2 = ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj);
   
   if(beta==1)
      return (x-xj) / sqrt( r2 + c*c );
   else   
      return   beta * (x-xj) * pow( r2 + c*c,   beta/2.0 - 1.0);   
}
//----------------------------------------------------------
template <typename T> 
T MQ<T>::dy(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2;
   
   r2 = ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj);
   
   if(beta==1)
       return (y-yj) / sqrt( r2 + c*c );
   else
       return  beta * (y-yj) * pow( r2 + c*c,  beta / 2.0 - 1.0);  
}
//----------------------------------------------------------
template <typename T> 
T MQ<T>::dz(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2;
   
   r2 = ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj);
   
   if(beta==1)
       return (y-yj) / sqrt( r2 + c*c );
   else
       return  beta * (z-zj) * pow( r2 + c*c,  beta / 2.0 - 1.0);  
}
//----------------------------------------------------------
template <typename T> 
T MQ<T>::dxx(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T  r2,c2;
   T  den1,den2;
   T  factor1;
   
   c2 = c*c;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj) + (z-zj)*(z-zj);
   
   if(beta==1){
        
        den1  = pow( r2 + c2 ,   3.0/2.0 );
        
        den2  = sqrt(r2+c*c);
        
        return  -((x-xj)*(x-xj))/den1 + 1.0/den2;
        
   }    
   else{
  	
  	   factor1 = pow( r2 + c2,  beta/2.0 - 1.0 ) ;
  	  
       return   beta * factor1 * ( 1.0 + ( ( beta-2.0) * ( (x-xj)*(x-xj) ) )/(r2 + c2) );    
   }    
}
//----------------------------------------------------------
template <typename T> 
T MQ<T>::dyy(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T  r2,c2;
   T  den1,den2;
   T  factor1;
   
   c2 = c*c;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj) + (z-zj)*(z-zj);
   
   if(beta==1){
        
        den1  = pow( r2 + c2 ,   3.0/2.0 );
        
        den2  = sqrt(r2+c*c);
        
        return  -((x-xj)*(x-xj))/den1 + 1.0/den2;
        
   }    
   else{
  	
  	   factor1 = pow( r2 + c2,  beta/2.0 - 1.0 ) ;
  	  
       return   beta * factor1 * ( 1.0 + ( ( beta-2.0) * ( (y-yj)*(y-yj) ) )/(r2 + c2) );    
   }
}
//----------------------------------------------------------
template <typename T> 
T MQ<T>::dzz(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T  r2,c2;
   T  den1,den2;
   T  factor1;
   
   c2 = c*c;
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj);
   
   if(beta==1){
        
        den1  = pow( r2 + c2 ,   3.0/2.0 );
        
        den2  = sqrt(r2+c*c);
        
        return  -((x-xj)*(x-xj))/den1 + 1.0/den2;
        
   }    
   else{
  	
  	   factor1 = pow( r2 + c2,  beta/2.0 - 1.0 ) ;
  	  
       return   beta * factor1 * ( 1.0 + ( ( beta-2.0) * ( (z-zj)*(z-zj) ) )/(r2 + c2) );    
   }
}


} // RBF namespace

#endif // _RBF_MQ_CPP_
