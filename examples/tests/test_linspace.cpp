/***************************************************************************
 *
 * Copyright (C) 2009   Jose Antonio Munoz Gomez
 *        
 * This file is part of Radial++
 * http://sourceforge.net/projects/radial/
 *
 * Radial++ is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 *
 * Radial++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For more information, please see the  Home Page:
 *    http://www.dci.dgsca.unam.mx/pderbf
 *
 ***************************************************************************/



#include "radial.h"



//---------------------------------------------------------
int main(int argc, char **argv)
{
  Vector<double>      x;
  
  x = linspace(10,10.,4);
  cout<<x<<endl<<endl;  
  
  x = linspace(0,0.,-2);
  cout<<x<<endl<<endl;
  
  x = linspace(1,0.,4);
  cout<<x<<endl<<endl;  
  
  x = linspace(0,1.,4);
  cout<<x<<endl<<endl;  
  
  x = linspace(1,0.,2);
  cout<<x<<endl<<endl;  
  
  return 0;
}
//---------------------------------------------------------





