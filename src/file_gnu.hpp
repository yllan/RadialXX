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
 
#ifndef  _FILE_GNU_HPP_
#define  _FILE_GNU_HPP_

namespace rbf{

#include <cstdio>

//----------------------------------------------------------
// Data 1D
//----------------------------------------------------------
template <typename T>
void save_gnu_data(const char *filename, T &x, T &y)
{
  FILE *  fpt = fopen(filename,"w");

   
 if(fpt==NULL)
 {
    fprintf(stderr,"ERROR: unabled to create file '%s'.\n\n",filename);
    exit(1);
 }

   
  for(int i=0; i < x.GetSize(); i++)
     fprintf(fpt,"%e  %e\n",x(i),y(i));     

fclose(fpt);
}   
//----------------------------------------------------------
// Data 2D
//----------------------------------------------------------
template <typename T>
void save_gnu_data(const char *filename, T &x, T &y, T &z)
{
  FILE  *  fpt = fopen(filename,"w");
   
   
 if(fpt==NULL)
 {
     fprintf(stderr,"ERROR: unabled to create file '%s'.\n\n",filename);
     exit(1);
 }   
   
  for(int i=0; i < x.GetSize(); i++)
   fprintf(fpt,"%e  %e  %e\n",x(i),y(i),z(i));     
   
 fclose(fpt);
}   
//----------------------------------------------------------
// Data 3D
//----------------------------------------------------------
template <typename T>
void save_gnu_data(const char *filename, T &x, T &y, T &z, T &f)
{
 FILE *  fpt = fopen(filename,"w");

   
  if(fpt==NULL)
  {
   fprintf(stderr,"ERROR: unabled to create file '%s'.\n\n",filename);
   exit(1);
  }   
   
  for(int i=0; i < x.GetSize(); i++)
   fprintf(fpt,"%e  %e  %e  %e\n",x(i),y(i),z(i),f(i));     
   
  fclose(fpt);
}
//----------------------------------------------------------
// Data n-dimensional
//----------------------------------------------------------
template <typename Vec>
void save_gnu_data(const char *filename, Vec &x, Vec &f, int dim)
{
 FILE *  fpt = fopen(filename,"w");

   
  if(fpt==NULL)
  {
   fprintf(stderr,"ERROR: unabled to create file '%s'.\n\n",filename);
   exit(1);
  }   
   
  for(int i=0; i < (x.GetSize()/dim); i++)
  {
    for(int j=0; j < dim; j++)
    {
         fprintf(fpt,"%e  ",x(i*dim + j));
    }
    fprintf(fpt,"%e\n",f(i));
  }
   
  fclose(fpt);
}   

} // RBF namespace

#endif //_FILE_GNU_HPP_




