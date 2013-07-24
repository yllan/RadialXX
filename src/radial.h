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
 
#ifndef _RADIAL_H_
#define _RADIAL_H_


#include <iostream>
#include <limits>
#include <cmath>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>


using namespace std;

//#define RADIAL_WITH_DEBUG   //show a naive information, only useful for the developer of Radial++
//#define RADIAL_WITH_C_CODE  //useful when no BLAS and LAPACK interface and low profile flags compilation


/*!
 *  Tool for time measurements (borrowed from FLENS).
 */
struct timer
{
/*!
 *  Start the clock.
 */
    void tic() { _time = clock(); }
/*!
 *  Stop the clock.
 */ 
    double toc() const { return double(clock() - _time)/double(CLOCKS_PER_SEC); }
    clock_t _time;
};


// Matrix and Vector classes from Seldon
#include "Seldon-5.0.1/Seldon.hxx"
using namespace Seldon;


// Polynomial required in the theory of radial basis functions
#include "pol/src/pol_polinomio.hpp"
using namespace pol;


// Radial Basis Functions
#include "rbf/rbf_all.h"
using namespace rbf;

//LU factorization for linear system of equations
#include "solver/solver_lu.hpp"
using namespace solver;


// Save files in ASCII format
#include "file_gnu.hpp"


// Make the uniform knot distribution 
#include "grid.hpp"


// Validate the fixed step-time in Euler method
#include "ode_control.hpp"


//The main routines to build the derivates matrices
//requiered in th unsymmetric collocation method
#include "fill/fill_select.hpp"
#include "fill/fillMatrix-1d.hpp"
#include "fill/fillMatrix-2d.hpp"
#include "fill/fillMatrix-3d.hpp"


//(x,y)   -> X=(x1,y1,x2,y2,...)
//(x,y,z) -> X=(x1,y1,z1,x2,y2,z2,...)
#include "compose.hpp"


//Data interpolation
#include "interpolation.hpp"


//Build the Gramm matriz for scattered data interpolation
#include "gram.hpp"


//Developed by Radial++ to Seldon
#include "radial2seldon/radial2seldon.h"


// to display data with  GNUPLOT
#ifdef    WITH_GNUPLOT
  #include "GNUplot.hpp" 
#endif


#endif //_RADIAL_H_
