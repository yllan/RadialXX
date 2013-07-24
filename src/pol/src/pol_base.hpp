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


#ifndef _POL_BASE_HPP_
#define _POL_BASE_HPP_


namespace pol{


template <typename T>
class Base
{
public:
  int     d;          // dimension del espacio
  T      constante; 
  int    *exponente;  //   X = constante * x_1 x_2 ... x_d  en R_d
                      //   exponentes de un punto en R_d
  bool initialized;
 
public:
   Base(void);
   Base(const Base &a);
  ~Base(void);
   void   free(void);
   void   make(int d);
   void   assign_values(T, int *exp,int d);
   void   deriva_parcial(int parcial);
   T      eval(const T *x, int d);
   Base&  operator=(const Base &p);
   void   show(void);
};

//---------------------------------------------------------------------
template<typename T>
Base<T>::Base(void)
{
             d = 0; 
   constante   = 0;
   exponente   = NULL;
   initialized = false;
}
//---------------------------------------------------------------------
template<typename T>
Base<T>::Base(const Base& a)
{
   int i;

          d   = a.d; 
  constante   = a.constante;
  initialized = a.initialized;
  
   exponente = new int[d];
   
   if(exponente==NULL){
       fprintf(stderr,"\nERROR: memory not created in the constructuor of Class Base in pol_base.cpp\n\n");
       exit(1);
    }
   
   for( i=0 ; i<d ; i++ )
      exponente[i] = a.exponente[i];
}
//---------------------------------------------------------------------
template<typename T>
void Base<T>::free(void)
{
   if(d>0){
      delete[] exponente;
   }  

   d = 0;
   initialized = false;
   exponente   = NULL;
}
//---------------------------------------------------------------------
template<typename T>
Base<T>::~Base(void)
{
   this->free();
}
//---------------------------------------------------------------------
 //recorro los elementos del polinomio derivado
template<typename T>
void Base<T>::make(int d)
{
   int  i;

   this->free();
 
   this->d   = d;
   exponente = new int[d];

   if(exponente==NULL){
       fprintf(stderr,"\nERROR: memory not created in Base.make(int).\n\n");
       exit(1);
   }

   initialized = true;
   constante   = 0;

   for( i=0; i<d; i++ )
   {
      exponente[i] = 0;
   }
}
//---------------------------------------------------------------------
template<typename T>
void Base<T>::assign_values(T constante, int *exp, int d)
{
   int i;
 
   if(initialized==false){
      fprintf(stderr,"\nWARNING: make not used after assign_values in the class Base.\n\n");
      make(d);
   }
   
   if(this->d!=d){
      fprintf(stderr,"\nERROR: dimension d=%d must be match with the used in function make.\n\n",d);
      exit(1);
   }

  this->constante = constante;

  for( i=0; i < d; i++)
  {
     exponente[i] = exp[i];
  }
}
//---------------------------------------------------------------------
template<typename T>
void Base<T>::show(void)
{
 if(d>0)
 {
  if(constante==0)
     return;
    
    printf("%.1f ",constante);
    
    switch(d){
       case 1: if(exponente[0]==0){
                }
               else
                printf("x^%d",exponente[0]);    
          break;
          
       case 2: if(exponente[0]==0){
               }
               else
                printf("x^%d",exponente[0]);    
            if(exponente[1]==0){
            }
            else
               printf("y^%d",exponente[1]);    
          break;
       case 3:if(exponente[0]==0){
               } 
             else if(exponente[0]==1)
               printf("x");
            else   
               printf("x^%d",exponente[0]);    
          
            if(exponente[1]==0){
             }
             else if(exponente[1]==1)
               printf("y");
              else
               printf("y^%d",exponente[1]);    
            if(exponente[2]==0){
            }
             else if(exponente[2]==1)
               printf("z");
             else
               printf("z^%d",exponente[2]);   
          break;
       default:
              for(int i=0;i<d;i++)
                printf("%d ",exponente[i]);
          
          
          break;
          
    }
    
    printf(", ");
 }
}
//---------------------------------------------------------------------
template<typename T>
void Base<T>::deriva_parcial(int ip)
{
   int exp;
 
//valido el elemento que quieren derivar
// X = constante * x_1 x_2 x_3 ... x_d
//
// ip: indica el subindice en el elemento X que quiero derivar
   if(ip<1){
                fprintf(stderr,"\nERROR: el indice del elemento a derivar debe ser al menos 1.\n\n");
                exit(1);
   }
  
   if(ip>d){
                fprintf(stderr,"\nERROR: el indice del elemento a derivar debe ser a lo mas %d.\n\n",d);
                exit(1);
   }
  
//valido si la constante es cero
   if(constante==0){
      for(int i=0;i<d;i++)
        exponente[i] = 0;
      
    return;   
   }
  
   ip = ip -1;
 // printf("exponente[ip] = %d , ip = %d\n",exponente[ip],ip);
//si el exponente de la variable que quiero derivar es uno, entonces
//al derivarlo todo vale cero
   if(exponente[ip]<=0){
      for(int i=0;i<d;i++)
        exponente[i] = 0;
      
    constante = 0;   
    return;   
   }
 
//derivo el elemento X en el indice ip
   exp           = exponente[ip]; 
   constante     = constante*exp;
   exponente[ip] = exp-1;

}
//----------------------------------------------------------------
template<typename T>
inline T pol_power(T x, int y)
{

        if(y==0)  return 1.0;
   else if(y==1)  return x;
   else if(y==2)  return x*x;
   else if(y==3)  return x*x*x;
   else if(y==4)  return x*x*x*x;
   else if(y==5)  return x*x*x*x*x; 
   else if(y==6)  return x*x*x*x*x*x; 
   else if(y==7)  return x*x*x*x*x*x*x; 
 
   return pow(x,(T)y);
}
//---------------------------------------------------------------------
template<typename T>
T Base<T>::eval(const T *x, int dimension)
{
   T aux,s,k;
   register int    i;
 
   if(this->d!=dimension){
     fprintf(stderr,"\nERROR: the dimension = '%d' in eval is different d=%d as the used in make.\n\n",dimension,this->d);
     exit(1);
   }

   if(constante==0)
   {
       return 0.0;
   }
  
   k    = constante;
   s    = 1.0;

   
   for( i=0; i<d; i++ )
   {
       aux = pol_power(x[i],exponente[i]);
       s   = s*aux;
   }

   return s*k;
}
//---------------------------------------------------------------------
//verbatim copy
template<typename T>
Base<T>& Base<T>:: operator=(const Base<T> &p)
{

   if(d!=p.d){
        delete[] exponente;

   exponente = new int[p.d];
   
      if(exponente==NULL){
        fprintf(stderr,"\nERROR: memory not created in the operator '=' in class Base.\n\n");
        exit(1);
      }
   }
  
              d = p.d; 
    constante    = p.constante;
    initialized  = p.initialized;
  
   
   for(int i=0;i<d;i++)
      exponente[i] = p.exponente[i];
    
    
   return *this;    
}   
//----------------------------------------------------------------
//No pertenece a la clase Polinomio
template<typename T>
Base<T> operator*(const T alpha, Base<T>  &p)
{
   return p.operator*(alpha);
}

} // namespace POL

#endif // _POL_BASE_HPP_


