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

#ifndef _MRBF_MQ_H_
#define _MRBF_MQ_H_

namespace rbf{

// beta = 1, 3, 5, 7, 9, 11, ...
// phi(r) =  ( 1 + c^2*r^2 )^(beta/2)   
template <typename T>
class MQM
{
     T   beta;
  bool   initialized;
  bool   initialized_pol;
   int    degree;
public:
  MQM(void);

  void         set_beta(int beta_factor);
  int          get_beta(void);   

//Determina el minimo grado del polinomio
  int          get_degree_pol(void);
  void         set_degree_pol(int degree_factor);
  int          get_min_degree_pol(void);   
   
  string       name(void){return "MQM";};      

//Evaluacion del MQM para datos en 1D, 2D, 3D
  T eval(T x,T xj, T c);
  T eval(T x,T y, T xj, T yj, T c);
  T eval(T x,T y, T z, T xj, T yj, T zj, T c);
   

  T eval(const T *x, const T *xj, int dim, T c);
};
   
} // RBF namespace

#endif // _RBF_MQ_H_
