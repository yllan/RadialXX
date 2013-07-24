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


#ifndef _RBF_TPS_H_
#define _RBF_TPS_H_

namespace rbf{
   
template <typename T>
class TPS
{
protected:   
    T    beta;
 bool    initialized;
 bool    initialized_pol;
  int    degree;   
public:
  TPS(void);

   void         set_beta(int beta_factor);
   int          get_beta(void);   
   
 //Determina el minimo grado del polinomio
   int          get_degree_pol(void);
   void         set_degree_pol(int degree_factor);
   int          get_min_degree_pol(void);   
   
   string       name(void){return "TPS";};   
   
//Evaluacion del nucleo radial para datos en 1D, 2D, 3D
  T eval(T x,T xj, T c);
  T eval(T x,T y, T xj, T yj, T c);
  T eval(T x,T y, T z, T xj, T yj, T zj, T c);

//new June - 2009
  T eval(const T *x, const T *xj, int dim, T c);   
  
  

//Derivatives
  T  dx(T x, T xj, T c);
  T  dxx(T x, T xj, T c);

  T  dx(T x, T y, T xj, T yj, T c);
  T  dy(T x, T y, T xj, T yj, T c);
  T  dxx(T x, T y, T xj, T yj, T c);
  T  dyy(T x, T y, T xj, T yj, T c);   
  T  dxy(T x, T y, T xj, T yj, T c);
  T  dyx(T x, T y, T xj, T yj, T c);  
  
  T  dx(T x, T y, T z, T xj, T yj, T zj, T c);
  T  dy(T x, T y, T z, T xj, T yj, T zj, T c);
  T  dz(T x, T y, T z, T xj, T yj, T zj, T c);
  T  dxx(T x, T y, T z, T xj, T yj, T zj, T c);
  T  dyy(T x, T y, T z, T xj, T yj, T zj, T c);  
  T  dzz(T x, T y, T z, T xj, T yj, T zj, T c); 
};

} // RBF namespace

#endif //_RBF_TPS_H_


