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
 /***************************************************************************
 * This file has been written to Radial++ in spirit to incorporate
 * the classical binary operations between vectors and matrices
 * in the library called Seldon.
 *
 * The original idea of this functions was taken from Blitz and suggested 
 * by Vivien Mallet.
 *
 ***************************************************************************/


#ifndef _MATH_FUNCTIONS_HPP_
#define _MATH_FUNCTIONS_HPP_


//
//! Vector <-- function ( Vector ) 
//!
//! Matrix <-- function ( Matrix ) 
//!
/*! 
    \param the 'functions' are definined below in this file 
*/
#define SELDON_DEFINE_FUNC(name)                     \
template <typename T>                                \
Vector<T> name (const Vector<T>& X)                  \
{                                                    \
   Vector<T> C( X );                                 \
   T *data   = C.GetData();                          \
   int m     = C.GetSize();                          \
                                                     \
   for( int i=0;  i<m;  i++ )                        \
      data[i] = name(data[i]);                       \
                                                     \
   return C;                                         \
}                                                    \
 template <typename T>                               \
Matrix<T> name (const Matrix<T>& X)                  \
{                                                    \
   Matrix<T> C( X );                                 \
   T *data   = C.GetData();                          \
   int m     = C.GetSize();                          \
                                                     \
   for( int i=0;  i<m;  i++ )                        \
      data[i] = name(data[i]);                       \
                                                     \
   return C;                                         \
}                                                    \



SELDON_DEFINE_FUNC(sin)
SELDON_DEFINE_FUNC(cos)
SELDON_DEFINE_FUNC(tan)
SELDON_DEFINE_FUNC(acos)
SELDON_DEFINE_FUNC(asin)
SELDON_DEFINE_FUNC(atan)

SELDON_DEFINE_FUNC(cosh)
SELDON_DEFINE_FUNC(sinh)
SELDON_DEFINE_FUNC(tanh)

SELDON_DEFINE_FUNC(exp)
SELDON_DEFINE_FUNC(log)
SELDON_DEFINE_FUNC(log10)


SELDON_DEFINE_FUNC(sqrt)

SELDON_DEFINE_FUNC(ceil)
SELDON_DEFINE_FUNC(fabs)
SELDON_DEFINE_FUNC(floor)


//Falta:  atan2 pow

//My own functions
template<typename T>
inline T pow2(const T& xc)
{
  return xc*xc;   
}
template<typename T>
inline T pow3(const T& xc)
{
  return xc*xc*xc;   
}
template<typename T>
inline T pow4(const T& xc)
{
  return xc*xc*xc*xc;   
}
template<typename T>
inline T pow5(const T& xc)
{
  return xc*xc*xc*xc*xc;   
}
template<typename T>
inline T pow6(const T& xc)
{
  return xc*xc*xc*xc*xc*xc;   
}
template<typename T>
inline T pow7(const T& xc)
{
  return xc*xc*xc*xc*xc*xc*xc;   
}
template<typename T>
inline T pow8(const T& xc)
{
  return xc*xc*xc*xc*xc*xc*xc*xc;   
}


SELDON_DEFINE_FUNC(pow2)
SELDON_DEFINE_FUNC(pow3)
SELDON_DEFINE_FUNC(pow4)
SELDON_DEFINE_FUNC(pow5)
SELDON_DEFINE_FUNC(pow6)
SELDON_DEFINE_FUNC(pow7)
SELDON_DEFINE_FUNC(pow8)




//
//! Vector <-- function ( Vector , Vector ) 
//!
/*! 
    \param the 'functions' are definined below in this file 
*/
#define SELDON_DEFINE_FUNC2(name)                    \
template <typename T>                                \
Vector<T> name (const Vector<T>& X,                  \
                const Vector<T>& Y)                  \
{                                                    \
   Vector<T> C( X.GetSize() );                       \
   int m     = C.GetSize();                          \
                                                     \
                                                     \
   CheckDims(X,Y,name);                              \
   for( int i=0;  i<m;  i++ )                        \
      data[i] = name(X(i),Y(i));                     \
                                                     \
   return C;                                         \
}                                                    \




//
//! Vector <-- Vector ^ alpha
//!
//! when alpha is a positive integer, it is more convenient to use the functions pow?, ?=2,3,..,8 
/*! 
    \param alpha scalar.
*/
template<typename T, typename Tc>
inline Vector<T> operator^(const Vector<T>& X, Tc alpha)
{
   Vector<T> C( X );
   T *data   = C.GetData();                          
   int m     = C.GetSize();                              

   for( int i=0;  i<m;  i++ )                        
      data[i]  = pow(data[i] , alpha );                   
      
   return C;
}    


//
//! Matrix <-- Matrix ^ alpha
//!
//! when alpha is a positive integer, it is more convenient to use the functions pow?, ?=2,3,..,8 
/*! 
    \param alpha scalar.
*/
template<typename T, typename Tc>
inline Matrix<T> operator^(const Matrix<T>& X, Tc alpha)
{
   Matrix<T> C( X );
   T *data   = C.GetData();                          
   int m     = C.GetSize();                              

   for( int i=0;  i<m;  i++ )                        
      data[i]  = pow(data[i] , alpha );                   
      
   return C;
}    

#endif //_MATH_FUNCTIONS_HPP_

