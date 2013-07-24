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


#ifndef _RBF_IMQ_CPP_
#define _RBF_IMQ_CPP_

#include "rbf_imq.h"

namespace rbf{
   
//----------------------------------------------------------
template <typename T> 
IMQ<T>::IMQ(void)
{
   beta            = 1; 
   degree          = get_min_degree_pol();
   initialized     = false;
   initialized_pol = false;
}
//----------------------------------------------------------------
template <class T>
std::ostream& operator<<(std::ostream &s,  IMQ<T> &rbf)
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
void IMQ<T>::set_beta(int beta_factor)
{
   //beta = 1, 3, 5, 7, 9, 11, ...
   // phi(r) =  (r^2 + c^2)^(-beta/2)
   
   if(beta_factor<1)
   {
      fprintf(stderr,"ERROR: in IMQ kernel, beta = %d debe ser un numero positivo impar 1,3,5,. \n\n",beta_factor);
      exit(1);
   }   
   
   
   if(!rbf_is_impar(beta_factor))
   {
      fprintf(stderr,"ERROR: beta = %d debe ser un numero impar. \n\n",beta_factor);
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
int IMQ<T>::get_beta(void)
{
  return int(beta);   
}
//----------------------------------------------------------
template<typename T>
int IMQ<T>::get_degree_pol(void)
{
  return degree;   
}
//----------------------------------------------------------
template <typename T>
void  IMQ<T>::set_degree_pol(int degree_factor)
{
  if(degree_factor<0){
    fprintf(stderr,"\n!! ERROR in IMQ-set_degree_pol: the degree = %d must be >= 0.\n\n",degree_factor);
    fflush(stderr);
    exit(1);
  }


   //First, validate that degre_factor >= min_degree_pol requiered
   if( degree_factor < get_min_degree_pol() )
   {
      fprintf(stderr,"WARNING: in IMQ kernel the degree = %d must be at least %d.\n\n",degree_factor, get_min_degree_pol());
   }   
   
   degree = degree_factor;   
   
   initialized_pol = true;   
}
//----------------------------------------------------------
template<typename T>
int IMQ<T>::get_min_degree_pol(void)
{
   return 0;   
}
//----------------------------------------------------------
template <typename T> 
T IMQ<T>::eval(T x,T xj, T c)
{
   T r;
   
   r = (x-xj)*(x-xj);      
   
 if(beta==1)
   return 1.0/sqrt( r + c*c ); 
 else
   return pow( r + c*c, -beta/2.0); 
  
}
//----------------------------------------------------------
template <typename T> 
T IMQ<T>::eval(T x, T y, T xj, T yj, T c)
{
   T r;
   
   r = (x-xj)*(x-xj) + (y-yj)*(y-yj);      
   
 if(beta==1)
   return 1.0/sqrt( r  + c*c); 
 else
   return pow( r  + c*c , -beta/2.0); 

}
//----------------------------------------------------------
template <typename T> 
T IMQ<T>::eval(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r;
   
   r = (x-xj)*(x-xj) + (y-yj)*(y-yj) +(z-zj)*(z-zj);   
   
   if(beta==1)
      return 1.0/sqrt( r + c*c); 
   else
      return pow(  r  + c*c , -beta/2.0 ); 
}
//----------------------------------------------------------
template <typename T> 
T IMQ<T>::eval(const T *x, const T *xj, int dim, T c)
{ 
   T r;
   
   r = rbf_norm_square(x,xj,dim);   
   
   if(beta==1)
      return 1.0/sqrt( r + c*c); 
   else
      return pow(  r  + c*c , -beta/2.0); 
}
//----------------------------------------------------------   
//              Data  1-D  
//----------------------------------------------------------
template <typename T>
inline T IMQ<T>::dx(T x, T xj, T c)
{  
  T r2;
   
  r2 = (x-xj) * (x-xj);
   
  return -beta * (x-xj) * pow( r2 + c*c, -beta/2.0 - 1.0  );
}
//----------------------------------------------------------
template <typename T>
inline T IMQ<T>::dxx(T x, T xj, T c)
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
inline T IMQ<T>::dx(T x, T y, T xj, T yj, T c)
{
   T r2;
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj); 
   
   return -beta * (x-xj) * pow( r2 + c*c, -beta/2.0 - 1.0  );   
   
}
//----------------------------------------------------------
template <typename T>
inline T IMQ<T>::dy(T x, T y, T xj, T yj, T c)
{
  T r2;
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj); 
   
   return -beta * (y-yj) * pow( r2 + c*c, -beta/2.0 - 1.0  );      
   
}
//----------------------------------------------------------
template <typename T>
inline T IMQ<T>::dxx(T x, T y, T xj, T yj, T c)
{
   T r2,c2;
   
   c2 = c*c;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj);
   
   return  -beta * pow( r2 + c2, -beta/2.0 - 1.0) * ( 1.0 + ( (-beta-2.0) * ( (x-xj)*(x-xj) ) )/(r2 + c2) );
   
}
//----------------------------------------------------------
template <typename T>
inline T IMQ<T>::dyy(T x, T y, T xj, T yj, T c)
{
   T r2,c2;
   
   c2 = c*c;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj);
   
   return  -beta * pow( r2 + c2, -beta/2.0 - 1.0) * ( 1.0 + ( (-beta-2.0) * ( (x-xj)*(x-xj) ) )/(r2 + c2) );   
}
//----------------------------------------------------------
template <typename T>
inline T IMQ<T>::dyx(T x, T y, T xj, T yj, T c)
{
   T  r2;   
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj) ;   
   
   return  -beta * (-beta-2.0) * (x-xj) * (y-yj) * pow( r2 + c*c,  -beta / 2.0 - 2.0);    	
}	
//----------------------------------------------------------
template <typename T>
inline T IMQ<T>::dxy(T x, T y, T xj, T yj, T c)
{
   T  r2;   
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj) ;   
   
   return  -beta * (-beta-2.0) * (x-xj) * (y-yj) * pow( r2 + c*c,  -beta / 2.0 - 2.0);   	
}	
//----------------------------------------------------------   
//              Data  3-D  
//----------------------------------------------------------   
template <typename T> 
T IMQ<T>::dx(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2;
   
   r2 = ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj);
   
   return   -beta * (x-xj) * pow( r2 + c*c,   -beta/2.0 - 1.0);   
}
//----------------------------------------------------------
template <typename T> 
T IMQ<T>::dy(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2;
   
   r2 = ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj);
   
   return  -beta * (y-yj) * pow( r2 + c*c,  -beta / 2.0 - 1.0);  
}
//----------------------------------------------------------
template <typename T> 
T IMQ<T>::dz(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T r2;
   
   r2 = ( x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj);
   
   return  -beta * (z-zj) * pow( r2 + c*c,  -beta / 2.0 - 1.0);  
}
//----------------------------------------------------------
template <typename T> 
T IMQ<T>::dxx(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T  r2,c2;
   T  factor1;
   
   c2 = c*c;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj) + (z-zj)*(z-zj);
   

  	   factor1 = pow( r2 + c2,  -beta/2.0 - 1.0 ) ;
  	  
       return   -beta * factor1 * ( 1.0 + ( ( -beta-2.0) * ( (x-xj)*(x-xj) ) )/(r2 + c2) );    

}
//----------------------------------------------------------
template <typename T> 
T IMQ<T>::dyy(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T  r2,c2;
   T  factor1;
   
   c2 = c*c;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj) + (z-zj)*(z-zj);
   

  	   factor1 = pow( r2 + c2,  -beta/2.0 - 1.0 ) ;
  	  
       return   -beta * factor1 * ( 1.0 + ( ( -beta-2.0) * ( (y-yj)*(y-yj) ) )/(r2 + c2) );    

}
//----------------------------------------------------------
template <typename T> 
T IMQ<T>::dzz(T x, T y, T z, T xj, T yj, T zj, T c)
{ 
   T  r2,c2;
   T  factor1;
   
   c2 = c*c;
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj) + (z-zj) * (z-zj);

  	
  	   factor1 = pow( r2 + c2,  -beta/2.0 - 1.0 ) ;
  	  
       return   -beta * factor1 * ( 1.0 + ( ( -beta-2.0) * ( (z-zj)*(z-zj) ) )/(r2 + c2) );    

}

} // RBF namespace

#endif // _RBF_IMQ_CPP_

