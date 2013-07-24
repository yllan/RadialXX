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

#ifndef _LOGSPACE_HPP_
#define _LOGSPACE_HPP_

namespace Seldon{

//----------------------------------------------------------
//
//     logspace(left, right , n )
//
//     Generate logarithmically spaced vectors
//     Similar to the implementation in Octave
//
//---------------------------------------------------------- 
template<typename T>
Vector<T> logspace(T xmin, T xmax, int n)
{
   Vector<T> x;   
   T   h;
   int i;
  
//if n<0  x=xmax
   if( n < 0 )
   {
         x.Reallocate(1);
         x(0) = pow( 10.0 , xmax );
         return x;
   }  
   
//if n==1   x=xmax
   if( n == 1 )
   {
         x.Reallocate(1);
         x(0) = pow( 10.0 , xmax );
         return x;      
   }  
   
//if n==2   x=[xmin , xmax]
   if( n == 2 )
   {
       x.Reallocate(2);
       x(0) = pow( 10.0 , xmin );
       x(1) = pow( 10.0 , xmax );
       return x;   
   } 
   
//if xmin == xmax   x.Fill(xmax)
   if (xmin == xmax )
   {
      x.Reallocate(n);
      x.Fill( pow( 10.0 , xmax ) );
      return x;   
   }   
   
//two options:  ascending order  xmin<xmax,  descinding order xmin>xmax

   h  = (xmax-xmin)/T(n-1);

//now fill with the data   
   x.Reallocate(n);
  
   x(0)   = pow( 10.0 , xmin );
   x(n-1) = pow( 10.0 , xmax );
   
   for( i=1; i<(n-1);  i++ )
   {
      x(i) = pow( 10.0 , xmin + h*i );
   }
   
   return x;
}


//----------------------------------------------------------
//
//     logspace(left, right , n )
//
//     Generate logarithmically spaced vectors
//     Similar to the implementation in Octave
//
//---------------------------------------------------------- 

template<typename T>
Vector<T> logspace(int xmin, T xmax, int n)
{
   return logspace(T(xmin), T(xmax), n );
}

template<typename T>
Vector<T> logspace(T xmin, int xmax, int n)
{
   return logspace(T(xmin), T(xmax), n );
}




} // namespace Seldon

#endif // _LOGSPACE_HPP_


