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

#ifndef _TRBF_VIRTUAL_H_
#define _TRBF_VIRTUAL_H_

namespace rbf{


template <typename T>
 class TRBF
{
 public:   
    string  version(void);
 
    virtual void         set_beta(int beta_factor)=0;
    virtual int          get_beta(void)=0;   

    virtual int          get_degree_pol(void)=0;
    virtual void         set_degree_pol(int degree_factor)=0;
    virtual int          get_min_degree_pol(void)=0;
   
    virtual string       name(void)=0;

    virtual T eval(T x,T xj, T c)=0;
    virtual T eval(T x,T y, T xj, T yj, T c)=0;
    virtual T eval(T x,T y, T z, T xj, T yj, T zj, T c)=0;

    virtual T eval(const T *x, const T *y, int dim, T c)=0;   
};

template<typename T>
string TRBF<T>::version(void)
{
 return "RBF version 0.42";
}

} // RBF namespace

#endif // _TRBF_VIRTUAL_H_


