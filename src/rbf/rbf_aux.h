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


#ifndef _RBF_AUX_H_
#define _RBF_AUX_H_

namespace rbf{
   
 int  rbf_maximo(int a, int b);
 int  rbf_minimo(int a, int b); 
 bool rbf_is_impar(unsigned int x);
 bool rbf_is_par(unsigned int x);
 
 template <typename T>  T rbf_pow(T r, int beta);
 template <typename T>  T rbf_norm_square(const T *x, const T *y, int n);
   
} // RBF namespace

#endif //_RBF_AUX_H_
