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

#ifndef _GRID_HPP_
#define _GRID_HPP_

namespace rbf{

//----------------------------------------------------------
// Data 1D for PDE
//---------------------------------------------------------- 
template<typename Vec>
void make_data(double xmin, double xmax, int nx,  Vec &x, int &ni, int &n)
{
  double hx;
 
   hx = (xmax-xmin)/double(nx-1);

   x.Reallocate(nx);
  
   ni = 0;
  
//interior nodes
   for(int i=1;i<(nx-1);i++)
   {
      x(ni) = xmin + hx*i;
      ni++;
   }

   n = ni;

//boundary nodes
   x(n)   = xmin;
   x(n+1) = xmax;
   
   n = n + 2;
}

//----------------------------------------------------------
// Data 2D for PDE in a box
//---------------------------------------------------------- 
template<typename Vec>
void make_data(double xmin, double xmax, double ymin, double ymax, int nx, int ny, 
               Vec &x, Vec &y, int &ni, int &n)
{
  double hx,hy;
 
   hx = (xmax-xmin)/double(nx-1);
   hy = (ymax-ymin)/double(ny-1);
  
   x.Reallocate(nx*ny);
   y.Reallocate(nx*ny);
  
   ni = 0;
  
//interior nodes  
   for(int j=1;j<(ny-1);j++)
   {
     for(int i=1;i<(nx-1);i++)
     {
       x(ni) = xmin + hx*i;
       y(ni) = ymin + hy*j;
        ni++;
     }
   }
  
    n = ni;
  
//boundary nodes
   for(int j=0;j<ny; j++)
   {
     for(int i=0;i<nx; i++)
     {
       if(i==0 || j==0 || i==(nx-1) || j==(ny-1) )
       {
         x(n) = xmin + hx*i;
         y(n) = ymin + hy*j;
          n++;
       } 
     }
   }  
  
}
//----------------------------------------------------------
// Data 2D for PDE in a circle centered at (0,0) with rad = 1.0
//----------------------------------------------------------  
template<typename T, typename Vec>
void make_data_circle(T radius, int nx, Vec &x, Vec &y, int &ni, int &n)
//  Nodes on the circle centered at (0,0)
//  the total number of nodes = nx * nx 
/*
   n = 7
      7x7           5x5     3x3   1x1
  o o o o o o o
  o o o o o o o  o o o o o
  o o o o o o o  o o o o o  o o o
  o o o o o o o  o o o o o  o o o  o
  o o o o o o o  o o o o o  o o o
  o o o o o o o  o o o o o
  o o o o o o o
      r = 1       r=0.66    r=0.33 r=0 
 
    inc_r = 2/(n-1) = 2/6 = 1/3
     
   n = 6
  o o o o o o
  o o o o o o  o o o o
  o o o o o o  o o o o  o o
  o o o o o o  o o o o  o o
  o o o o o o  o o o o
  o o o o o o
   r = 1         r=0.66  r=0.33

   inc_r = 2/n = 1/3
    n = 8
   o o o o o o o o
   o             o  o o o o o o
   o             o  o         o  o o o o
   o             o  o         o  o     o  o o 
   o             o  o         o  o     o  o o
   o             o  o         o  o o o o
   o             o  o o o o o o
   o o o o o o o o
      r=1             r=0.75      r=0.5   r=0.25
      
      inc_r = 2/n = 2/8 = 1/4
      
     n = 10 
   o o o o o o o o o o
   o                 o
   o                 o
   o                 o
   o                 o
   o                 o
   o                 o   8x8     6x6    4x4    2x2   1*1
   o                 o
   o                 o    
   o o o o o o o o o o
     r = 1               r=0,8    r=0.6  r=0.4 r=0.2  r=0.0
     inc_r = 2/n = 2/10 = 1/5  
*/  
{
   int    m;
   int    n_nodes;
   double h;
   double xc,yc;
   double t,inc_r;
   double r;
   int    nx_original;
   int    cont=0;

   
   nx_original = nx;
   
//check the radius
   if(radius<1e-5){
      printf("ERROR: the radius is too small %f\n",radius);
      printf("       in file '%s' near to %d\n\n",__FILE__,__LINE__);
      exit(1);      
   }
   
//check the minimum number of nodes
  if(nx<=3){
       printf("WARNING: the grid %d x %d is too small.\n\n",nx,nx);
  }   


   r = radius;

// check if n is impar
   if((nx%2)==0)
      inc_r = 2.0/nx;
   else   
      inc_r = 2.0/(nx-1);
      
//determine the total number of nodes
   if((nx%2)==0)
       m = nx*nx + 1;
   else
       m = nx*nx;     
      
//resize the output vectors   
   x.Reallocate(m);
   y.Reallocate(m);
   
//external loop for interior nodes
   while(nx>1)
   {
      
      n_nodes =  nx + nx + nx-2 + nx-2;
      h       = (2.0*M_PI)/(n_nodes);
      
      //printf("n_nodes = %d\n",n_nodes);
      //printf("nx      = %d\n",nx);
      //printf("radius  = %f\n",r);
      
      for(int i=0; i<n_nodes; i++)
      {
              t    = 0.0 + i*h;
            xc   = r*cos(t);
            yc   = r*sin(t);   
         if(r<radius)
         {   
            x(cont) = xc;
            y(cont) = yc;
            cont++;
            
            //printf("%f %f\n",xc,yc);
         }  
      }   
        
     nx = nx - 2;   
     r  = r - inc_r;
   }   

//add the center of the circle
   x(cont) = 0.0;
   y(cont) = 0.0;
    cont++;
  
//assign the number of interior nodes  
   ni = cont; 
   
//external loop for boundary nodes   
   nx = nx_original;
   r  = radius;
      n_nodes =  nx + nx + nx-2 + nx-2;
      h       = (2.0*M_PI)/(n_nodes);
         
      for(int i=0; i<n_nodes; i++)
      {
              t    = 0.0 + i*h;
            xc   = r*cos(t);
            yc   = r*sin(t);   
            x(cont) = xc;
            y(cont) = yc;
            cont++;
      }   
      
//assign the total number of nodes: interior nodes + boundary nodes      
   n = x.GetSize();

}

} // RBF namespace

#endif // _GRID_HPP_



