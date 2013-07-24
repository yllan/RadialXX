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

#ifndef _RBF_MQ_DERIVATIVES_HPP_
#define _RBF_MQ_DERIVATIVES_HPP_

namespace rbf{


//----------------------------------------------------------   
//              Data  1-D  
//----------------------------------------------------------
template <typename T>
inline T mq_1d_x(int beta , T x, T xj, T c)
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
inline T mq_1d_xx(int beta , T x, T xj, T c)
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
inline T mq_2d_x(int beta , T x, T y, T xj, T yj, T c)
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
inline T mq_2d_y(int beta , T x, T y, T xj, T yj, T c)
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
inline T mq_2d_xx(int beta , T x, T y, T xj, T yj, T c)
{
   T  r2,c2;
   T  den1,den2;
   T  factor1;
   
   c2 = c*c;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj);
   
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
T mq_2d_yy(int beta, T x, T y, T xj, T yj, T c)
{
   T r2,c2;
   T den1,den2;
   
   c2 = c*c;
   
   r2 = (x-xj)*(x-xj) + (y-yj)*(y-yj);
   
   if(beta==1){
        
        den1   = pow( r2 + c2 ,  3.0/2.0 );
        
        den2   = sqrt(r2+c*c);
        
        return  -((y-yj)*(y-yj))/den1 + 1.0/den2;
        
   }   
   else   
      return   beta * pow( r2 + c2,  beta/2.0 - 1.0 ) * ( 1.0 + ( (  beta - 2.0) * ( (y-yj)*(y-yj) ) )/(r2 + c2) );   
   
}
//----------------------------------------------------------
template <typename T>
inline T mq_2d_yx(int beta, T x, T y, T xj, T yj, T c)
{
   T  r2;   
   
   r2 = (x-xj) * (x-xj) + (y-yj) * (y-yj) ;   
   
   return  beta * (beta-2.0) * (x-xj) * (y-yj) * pow( r2 + c*c,  beta / 2.0 - 2.0);    
}
//----------------------------------------------------------
template <typename T>
inline T mq_2d_xy(int beta, T x, T y, T xj, T yj, T c)
{
   return mq_2d_yx(beta, x, xj, y, yj, c);
}
//----------------------------------------------------------   
//              Data  3-D  
//----------------------------------------------------------   
template <typename T> 
T mq_3d_x(int beta , T x, T y, T xj, T yj, T z, T zj, T c)
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
T mq_3d_y(int beta , T x, T y, T xj, T yj, T z, T zj, T c)
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
T mq_3d_z(int beta , T x, T y, T xj, T yj, T z, T zj, T c)
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
T mq_3d_xx(int beta , T x, T y, T xj, T yj, T z, T zj, T c)
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
T mq_3d_yy(int beta , T x, T y, T xj, T yj, T z, T zj, T c)
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
T mq_3d_zz(int beta , T x, T y, T xj, T yj, T z, T zj, T c)
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

#endif // _RBF_MQ_DERIVATIVES_HPP_

