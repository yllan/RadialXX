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


#ifndef _POL_POLINOMIO_CPP_
#define _POL_POLINOMIO_CPP_


#include "pol_base.hpp"
#include "pol_bitvector.hpp"



namespace pol{

template<typename T> 
class Polinomio
{
protected:
    Base<T>       *data;
   int            d,m,M;
   bool     initialized;
   
public:
   Polinomio(void);   
  ~Polinomio(void);
   void    free(void);
   void    make(int d,int m);
   void    deriva(int ip, int orden);
   void    deriva(string ip, int orden);

   void      eval(const T *x, int d, T *px, int dim_px);
   T**       build(const T *x, int n, int d);
   Matrix<T> build_tnt(const T *x, int n, int d);
    
   Polinomio&  operator=(const Polinomio &a);
   
   int    get_M(void);
   int    get_d(void);
   int    get_m(void);
   void   show(void);
   
private:
   double factorial(double n);   
   int    combinaciones(int r, int c);
   void   deriva_parcial(int ip);
   void   make_data_1d(void);   
   void   make_data_2d(void);   
   void   make_data_3d(void);  
   void   make_data_nd(void);  
};
//----------------------------------------------------------------
template<typename T>
Polinomio<T>::Polinomio(void)
{
        d=m=M = 0;
  initialized = false;
         data = NULL;

}
//----------------------------------------------------------------
template<typename T>
void Polinomio<T>::free(void)
{
  if(initialized)
  {
     delete[] data;
  } 
   
        d=m=M = 0;
  initialized = false;
  data        = NULL;   
}
//----------------------------------------------------------------
template<typename T>
Polinomio<T>::~Polinomio(void)
{
   this->free();	
}
//----------------------------------------------------------------
template<typename T>
double Polinomio<T>::factorial(double n)
{
 int i;
 double    fact = 1;
  
   for(i = 1; i<= int(n); i++)
        fact = fact * (double)i;
    
 return fact;
}
//----------------------------------------------------------------
template<typename T>
int Polinomio<T>::combinaciones(int r, int c)
{
   return int(factorial(double(r))/(factorial(double(c))*factorial(double(r-c))));
}
//----------------------------------------------------------------
template<typename T>
void Polinomio<T>::make_data_1d(void)
{
   int  cont = 0;
   int  i;
   int  exp[3];
   
//Inicializo los elementos constantes de cada elemento de la base
   for( i = 0; i < M; i++ ) 
   {
      exp[0] = i;
               
      data[i].assign_values( 1.0 , exp, this->d);

      cont++;

      if(cont>M){
           fprintf(stderr,"\nERROR: the base elements %d exceed the calculated M = %d\n",cont,M);
           fprintf(stderr,"        in Polinomio.make 1d\n\n");
           exit(1);
      }
   }
   
   if(cont<M){
      fprintf(stderr,"\nERROR: the base elements %d are not equal than  M = %d\n",cont,M);
      fprintf(stderr,"        in Polinomio.make 1d\n\n");
      exit(1);   	
   }  
}      
//----------------------------------------------------------------
template<typename T>
void Polinomio<T>::make_data_2d(void)
{
   int  cont;
   int  i,j;
   int  exp[3];
   
   cont = 0;
//Inicializo los elementos constantes de cada elemento de la base
   for( i = 0; i < m; i++ )
   { 
      for( j = 0; j < m-i; j++ ) 
      {
         exp[0] = i;
         exp[1] = j;
               
         data[cont].assign_values( 1.0 , exp, this->d);

         cont++;

         if(cont>M){
           fprintf(stderr,"\nERROR: the base elements %d exceed the calculated M = %d\n",cont,M);
           fprintf(stderr,"        in Polinomio.make 2d\n\n");
           exit(1);
         }
      }
   }      
   
   if(cont<M){
      fprintf(stderr,"\nERROR: the base elements %d are not equal than  M = %d\n",cont,M);
      fprintf(stderr,"        in Polinomio.make 2d\n\n");
      exit(1);   	
   }   
           
}
//----------------------------------------------------------------
template<typename T>
void Polinomio<T>::make_data_3d(void)
{
   int  cont = 0;
   int  i,j,k;
   int  exp[3];


   for( k=0; k<m;k++)
      for( i=0; i<m-k;i++)
         for(j=0; j<m-i-k;j++)
         {
            exp[0] = k;
            exp[1] = i;
            exp[2] = j;
               
            data[cont].assign_values(1.0,exp,this->d);
            cont++;
            
            if(cont>M){
                fprintf(stderr,"\nERROR: the base elements %d exceed the calculated M = %d\n",cont,M);
                fprintf(stderr,"        in Polinomio.make 3d\n\n");
                exit(1);
            }            
            
         }
         
   if(cont<M){
      fprintf(stderr,"\nERROR: the base elements %d are not equal than  M = %d\n",cont,M);
      fprintf(stderr,"        in Polinomio.make 3d\n\n");
      exit(1);   	
   }  
}         
//----------------------------------------------------------------
template<typename T>
void Polinomio<T>::make_data_nd(void)
{
  TVectorBits vbits;
  int           sum;   
  int          *vec;   
  int          cont;
  
   cont = 0;

//create the dynamic memory
   vec = new int[d]; 
   
   if(vec==NULL){
      fprintf(stderr,"\nERROR: memory not allowed in Polinomio-make_nd.\n\n");
      exit(1);
   }
       
//create the base vector of longitud d in the base number m
   vbits.initialize(d,m);
            
//make the firsth entrance                       
   vbits.first(vec);
        
    
   do{
       sum = 0;
        
       //no problem with the sum, constains integers not to big     
         for(int i=0; i<d; i++)  
             sum  +=  vec[i];
            
        if(sum<m){
            
               data[cont].assign_values(1.0,vec,d);
               
               cont=cont+1;
               
               
               if(cont>M){
                    fprintf(stderr,"\nERROR: the base elements %d exceed the calculated M = %d\n",cont,M);
                    fprintf(stderr,"        in Polinomio.make nd\n\n");
                    exit(1);
               }      
               
         }
               
    }while(vbits.next(vec));
    
   delete[] vec;
   
   if(cont<M){
      fprintf(stderr,"\nERROR: the base elements %d are not equal than  M = %d\n",cont,M);
      fprintf(stderr,"        in Polinomio.make nd\n\n");
      exit(1);   	
   }  
}             
//----------------------------------------------------------------
template<typename T>
void Polinomio<T>::make(int d,int m)
{

  if(m<0){
    fprintf(stderr,"\nERROR in Polynomial-make: the degree m=%d must be >= 0.\n\n",m);
    exit(1);
  }
  
  this->free();
  
  if(m==0){   //m indica el orden del polinomio
              //indicate null polynomial
     
     initialized = true;
         this->d = d; // dimension
         this->M = 0; //numero de elementos de la base del polinomio
         this->m = m; //degree of the polynomial
         this->data = NULL;
     return;
  }

  if(m>10){
    fprintf(stderr,"\nINTERNAL-ERROR: too big the degree of the polynomial, degree = %d.\n",m);
    exit(1);
  }
  
  this->d    = d;
  this->m    = m; 
  this->M    = combinaciones(d+m-1,m-1);


//Reservo la memoria para el arreglo de elementos del polinomio
   data = new Base<T>[this->M];
   
   if(data==NULL){
      fprintf(stderr,"\nERROR: memory not allowed in Polinomio-make.\n\n");
      exit(1);
   }
   
//For each base element in pol, we reserve the memory, recall p(x)={p1,p2,...,pm}
//where pi is of type pol_base.hpp 
   for(int i=0; i<this->M; i++)
   {
      data[i].make(d);
   } 


   switch(d)
   {
     case 1:    make_data_1d();
         break;
          
     case 2:    make_data_2d();
         break;

     case 3:    make_data_3d();
         break;
         
     default:   make_data_nd();
         break;
   }
   
   initialized=true;
}
//---------------------------------------------------------------------
template<typename T>
void Polinomio<T>::deriva_parcial(int ip)
{
   for(int i=0; i<M; i++)
     data[i].deriva_parcial(ip);
}
//---------------------------------------------------------------------
template<typename T>
void   Polinomio<T>::deriva(int ip, int orden)
{
   for(int i=0; i<orden;i++)
     deriva_parcial(ip);
}
//---------------------------------------------------------------------   
template<typename T>
void  Polinomio<T>::deriva(string ip, int orden)
{
   if(ip=="x")
   {
      for(int i=0; i<orden;i++)
         deriva_parcial(1);
   }
   else if(ip=="y")
   {
      for(int i=0; i<orden;i++)
         deriva_parcial(2);
   }
   else if(ip=="z")
   {
      for(int i=0; i<orden;i++)
         deriva_parcial(3);
   }
   else{
     fprintf(stderr,"\nERROR: no esta definida la derivacion del polinomio con respecto de: '%s'.\n\n",ip.c_str());
     exit(1);
   }
}
//----------------------------------------------------------------
template<typename T>
void Polinomio<T>::eval(const T *x, int d, T *px, int dim_px)
{
	
   if(this->initialized==false){
     fprintf(stderr,"\nERROR: the polynomial has not been created, in 'eval'.\n\n");
     exit(1);
   }

   if( this->M!=dim_px ){
     fprintf(stderr,"\nERROR: en 'eval'  M='%d' debe coincidir con el argumento de entrada dim_lambda='%d'.\n\n",this->M,dim_px);
     exit(1);
   }

//evaluo el polinomio y cada elemento de la base del polinomio 
   for( int i=0; i<M; i++)
   {
      px[i] = data[i].eval(x,d);
   }
}
//----------------------------------------------------------------
template<typename T>
T** Polinomio<T>::build(const T *x, int n, int d)
// RowMajor   
{
   T **P;
   
   if(this->initialized==false){
      fprintf(stderr,"\nERROR: the polynomial has not been created, in 'build'.\n\n");
      exit(1);
   }

   if(this->M == 0)
      return NULL;

   P = new T*[n];


   for(int i=0; i<n; i++)
   {
       P[i] = new T[this->M];
   }

   for(int i=0; i<n; i++ )
   {
       eval( &x[i*d], d, P[i], M);
   }

  return P;
}
//---------------------------------------------------------------------   
template<typename T>   
Matrix<T>  Polinomio<T>::build_tnt(const T *x, int n, int d)
// RowMajor   
{
   Matrix<T>   P;   

   if(this->initialized==false)
   {
         fprintf(stderr,"\nERROR: the polynomial has not been created, in 'build'.\n\n");
         exit(1);
   }
      
   if(this->M == 0)
      return P;
      
   P.Reallocate(n,M);   
      
   T *y = new T[M];
      
   for(int i=0; i<n; i++ )
   {
       eval( &x[i*d], d, y, M);
         
       for( int j=0;  j<M;  j++ )
         P(i,j) = y[j];    
   }
   
   delete[] y;
      
  return P;
}   
//---------------------------------------------------------------------
//verbatim copy
template<typename T>
Polinomio<T>&  Polinomio<T>::operator=(const Polinomio<T> &p)
{

//checar su p ya se inicializado con make
   if(p.initialized==false){
      fprintf(stderr,"\nERROR: the polynomial has not been created, in operator '='.\n\n");
      exit(1);
   }   
   
   
   if(this->M!=p.M ||  this->m != p.m || this->d != p.d ){
       delete[]  data; // pol es del tipo PonyPol[]
       data  = new Base<T>[p.M];
       
   
       for(int i=0; i<p.M; i++)
          data[i].make(p.d);
      
       
       
       d    = p.d;
       M    = p.M;
       m    = p.m;

       if(data==NULL){
         fprintf(stderr,"\nERROR: memory not allowed in the function '=' in Polynomial Class.\n\n");
         exit(1);
       }
   }

   for(int i=0;i<M;i++)
      data[i] = p.data[i];


   this->d = p.d;
   this->m = p.m;
   this->M = p.M;   
   this->initialized = p.initialized;
   
   return *this;    
}
//---------------------------------------------------------------------
template<typename T>
void   Polinomio<T>::show(void)
{
   if(M>0)
   {
      printf("---------------------------------\n");
   
      for(int i=0; i<M; i++)
          data[i].show();
   
      printf("\n");  
   }
}
//----------------------------------------------------------------
template<typename T>
int Polinomio<T>::get_M(void)
{
   if(this->initialized==false){
     fprintf(stderr,"\nERROR: the polynomial has not been created, in 'get_M'.\n\n");
     exit(1);  
   }
  return M;
}
//----------------------------------------------------------------
template<typename T>
int Polinomio<T>::get_d(void)
{
   if(this->initialized==false){
      fprintf(stderr,"\nERROR: the polynomial has not been created, in 'get_d'.\n\n");
      exit(1);  
   }
   
  return d;
}
//----------------------------------------------------------------
template<typename T>
int Polinomio<T>::get_m(void)
{

   if(this->initialized==false){
   	 fprintf(stderr,"\nERROR: the polynomial has not been created, in 'get_m'.\n\n");
     exit(1);  
   }
   
  return m;
}
//----------------------------------------------------------------
template <class T>
std::ostream& operator<<(std::ostream &s,  Polinomio<T> &pol)
{
    s <<"Pol info: "<<"\n";
    s <<" data dimension : " << pol.get_d()<< "\n";
    s <<" order          : " << pol.get_m()<< "\n";
    s <<" base elements  : " << pol.get_M()<< "\n"; 
    return s;
}
//----------------------------------------------------------------
template <class T>
std::ostream& operator<<(std::ostream &s,  const Polinomio<T> &pol)
{
    s <<"Pol info: "<<"\n";
    s <<" data dimension : " << pol.get_d()<< "\n";
    s <<" order          : " << pol.get_m()<< "\n";
    s <<" base elements  : " << pol.get_M()<< "\n"; 
    return s;
}
   


} // namespace POL

#endif // _POL_POLINOMIO_CPP_


