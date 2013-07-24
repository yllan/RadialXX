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
 
//
// Interpolate data n-dimensional
//
// usage:
//         ./interpolate  example_1d.cfg
//

#ifdef WITH_LAPACK  
  #define SELDON_WITH_CBLAS
  #define SELDON_WITH_LAPACK
#endif  
#include "radial.h"

#define SIZE_NAME 512

//---------------------------------------------------------
int get_dim(char *filename)
// Taken from Seldon
{
   Vector<double> first_row;
   ifstream FileStream;
   string   line;
   
//open the file      
   FileStream.open(filename);
  
//read the first line.
   getline(FileStream, line);

//check if was opened
   if (FileStream.fail())
      return 0;

//converts the first line into a vector.
   istringstream line_stream(line);
 
   first_row.ReadText(line_stream);
    
//close the file
   FileStream.close();
   
//return the size of the vector   
   return first_row.GetSize();
}    
//---------------------------------------------------------
template <class RBF>
void interpolate(RBF rbf           ,
              int  dim             , 
              int  beta_band       , int beta_factor,
              int  pol_degree_band , int pol_degree,
              int  c_value_band    , double c_value,
              char* input_centers  ,
              char* input_data     ,
              char* output_data)
{
   Vector<double>    x,b,X,x_new,f;
   int               n,m;
   int               cont,cont_b;
   Polinomio<double> pol;
   Matrix<double>    A;
   Vector<double>    lambda;
   LU<double>        lu;
   
//optional: set the beta factor in thr rbf
   if(beta_band==1)
      rbf.set_beta(beta_factor);
    
//optional: set the exponent in the RBF 
   if(pol_degree_band==1) 
      rbf.set_degree_pol(pol_degree);
     
//configure the associate polynomial
// pol.make( data_dimension, degree_pol)  
   pol.make(dim , rbf.get_degree_pol());   
   
//show rbf and pol info
   cout<<rbf;
   cout<<pol<<endl;   
     
//load in memory the input_centers file
   X.ReadText(input_centers);
   
//load in memory the input_data file    
   x_new.ReadText(input_data); 
   
//obtain the number of center nodes that will be used to build the Gramm matrix   
    n = X.GetSize() /( dim +1);   
     
//validate that each file contain some data
   if(X.GetSize()==0)
   {
      fprintf(stderr,"ERROR: we not found elements in file '%s'.\n",input_centers);
      fprintf(stderr,"\n");
      exit(1);      
   } 
   if(x_new.GetSize()==0)
   {
      fprintf(stderr,"ERROR: we not found elements in file '%s'.\n",input_data);
      fprintf(stderr,"\n");
      exit(1);      
   } 
   
//validate the size of the vector and the data dimension
   if((X.GetSize()%(dim+1))!=0) 
   {
      fprintf(stderr,"ERROR: the elements size in '%s' are not divisible by dim + 1 = %d\n",input_centers,dim+1);
      fprintf(stderr,"       check that all the rows have the same elements.\n\n");
      exit(1);
   }
   
//validate the size of the vector and the data dimension   
   if((x_new.GetSize()%dim)!=0) 
   {
      fprintf(stderr,"ERROR: the elements size in '%s' are not divisible by dim = %d\n\n",input_data,dim);   
      fprintf(stderr,"       check that all the rows have the same elements.\n\n");
      exit(1);
   }
   
//validate the data-dimension   
   if( (get_dim(input_centers)-1) != dim)
   {
      fprintf(stderr,"ERROR: the data dimension detected in '%s' was = %d \n",input_centers,get_dim(input_centers)-1);
      fprintf(stderr,"       the data_dimension as input argument are = %d\n",dim);   
      fprintf(stderr,"       check the data dimension in the input configuration file\n");
      fprintf(stderr,"\n");
      exit(1);      
   }
   
//validate the data-dimension  
   if( (get_dim(input_data)) != dim)
   {
      fprintf(stderr,"ERROR: the data dimension detected in '%s' was = %d \n",input_centers,get_dim(input_centers)-1);
      fprintf(stderr,"       the data_dimension as input argument are = %d\n",dim);   
      fprintf(stderr,"       check the data dimension in the input configuration file\n");
      fprintf(stderr,"\n");
      exit(1);        
   }
      
//show info of the readed files
   printf("Number of centers = %d with dimension = %d in '%s'.\n",n,dim,input_centers);
   printf("\n");      
   printf("Number of nodes to interpolate = %d with dimension = %d in '%s'.\n",x_new.GetSize()/dim,dim,input_data);
   printf("\n");      
  
  
//validate that the number of nodes n in file 'input_centers' must be al least m-nodes
//the value of m correspond to the number of elements in the polynomial base.
   if(n<pol.get_M())
   {
      fprintf(stderr,"ERROR: the number of nodes n = %d in '%s' must be at least %d\n",n,input_centers,pol.get_M());
      fprintf(stderr,"       the value %d correspond to the number of elements in the base of the polynomial choosed\n",pol.get_M());   
      fprintf(stderr,"\n");
      exit(1);              
   } 
   
//extract from X the center nodes and the valuez of z=f(x1,..,x_dim)
   cont = 0; cont_b=0;
   x.Reallocate( n*dim );
   b.Reallocate(n);
   for(int i=0; i<X.GetSize(); i += dim + 1)
   {
       for(int j=0; j<dim; j++)
       {
         x(cont++) = X( i + j );   
       }
       b(cont_b++) = X( i+dim );
   }   
   
//clean the X vector, we liberate the memory
   X.Reallocate(0);
   

//get the number of elements in the polynomial base  
   m = pol.get_M();
   
   
//reallocate b without lost of the previously stored data   
   b.Resize(n+m);
   
 
//fill with zeros the region for the polynomial
//  A  P  = b
// P^t 0    0
   for( int i=0;  i<m;  i++)
       b(n+i) = 0.0;
   
 
//fill the Gramm matrix    with the rbf and the polynomial evaluations 
   printf("Building the Gramm matrix\n");
   printf("\n");
   fill_gram(rbf ,  pol, c_value, x, dim, A);   
  
//solve the linear system by LU factorization
   cout<<"Running the LU matrix factorization ( "<<A.GetM()<<" x "<<A.GetN()<<" ) "<<endl;
   printf("\n");
   lambda = lu.gauss(A, b); 
   
//interpolate the data   
   Vector<double> myc(n);
   
   myc.Fill( c_value );
      
   printf("\n");         
   printf("Running the interpolation process ( n = %d , m = %d )\n",n,x_new.GetSize()/dim);
   printf("\n");   
   
   f = interpolate(rbf,pol, myc , lambda,  x,x_new,dim); 

//save to file the interpolated data
   printf("Storing the interpolated data in '%s'\n",output_data);
   printf("\n");
   save_gnu_data(output_data,x_new,f,dim);

//end of the interpolation
   printf("\n");
   printf("End of the interpolation with Radial++ \n");
   printf("       J. Antonio Munoz-Gomez          \n");
   printf("\n");
}                   
                                      
//---------------------------------------------------------
int process_line(char* data)
// return 1 if data[0]=='#'
//        0 otherwise  
{
   
   if(strlen(data)<=1)
     return 1;   
    
   if(data[0]=='#')
     return 1;
     
  return 0;      
}
//---------------------------------------------------------
void remove_spaces_enter(char *token)
{
 if(token==NULL)
    return;
    
 int len = strlen(token);
 
 if(token[len-1]=='\n')
    token[len-1]='\0';       
}
//---------------------------------------------------------
void read_token(const char *filename,const char* token, char* token_output)
{
   FILE *in;
   char  data[512];
   char  *ptr;
   char  *tok1,*tok2,*tok3;
   
   
   in = fopen(filename,"r");
   
   if(in == NULL){
     fprintf(stderr,"\nERROR: unabled to open the filename '%s'.\n\n",filename);
     exit(1);   
   }
   
   strcpy(token_output,"");
   
   while(!feof(in))
   {
       ptr = fgets(data,511,in);
       
       if(ptr!=NULL){
          //now, process the input line
            if(process_line(data)==0)
            {
               //get three tokens
                  tok1 = strtok(data," ");
                  tok2 = strtok(NULL," ");
                  tok3 = strtok(NULL," ");    
                  
                //remove the posible enter '\n'
                  remove_spaces_enter(tok1);
                  remove_spaces_enter(tok2);
                  remove_spaces_enter(tok3);
               
               // check the token
                  if(strcmp(tok1,token)==0)
                  {
                     strcpy(token_output,tok3); 
                  }  
            }
       }
   }  
   
  fclose(in); 
  
}


//---------------------------------------------------------
int main(int argc, char **argv)
{
   FILE   *in;
   char   filename[SIZE_NAME];
   char   rbf[SIZE_NAME];
   char   beta[SIZE_NAME];
   char   shape_parameter[SIZE_NAME];
   char   polynomial_degree[SIZE_NAME];
   char   input_centers[SIZE_NAME];
   char   data_dimension[SIZE_NAME];
   char   input_data[SIZE_NAME];
   char   output_data[SIZE_NAME];
   int    rbf_type;
   double c_value;
   int    pol_degree;
   int    beta_factor;
   int    beta_band;
   int    pol_degree_band;
   int    c_value_band;
   int    dim;

//check the number of input arguments 
   if(argc<2){
     printf("%s  <filename>\n",argv[0]);
     printf("                example: %s example_1d.cfg\n",argv[0]);
     printf("                example: %s example_2d.cfg\n",argv[0]);
     printf("\n");
     return 1;
   }
   
//show info
   printf("\n");   
   printf(" Scattered Data Interpolation with Radial Basis Functions\n");
   printf("                    Radial++                             \n");
   printf("\n");   
   
//assign the filename with the configuration
   strcpy(filename,argv[1]);
   
//read the tokens
   read_token(filename,"RBF",rbf); 
   read_token(filename,"BETA",beta);    
   read_token(filename,"POLYNOMIAL_DEGREE",polynomial_degree);      
   read_token(filename,"SHAPE_PARAMETER",shape_parameter);
   read_token(filename,"FILE_INPUT_CENTERS",input_centers);
   read_token(filename,"FILE_INPUT_DATA",input_data);
   read_token(filename,"FILE_OUTPUT_DATA",output_data);
   read_token(filename,"DATA_DIMENSION",data_dimension);
   
//show the configuration readed   
   if(strlen(rbf)>0)               printf("  RBF                <--  %s\n",rbf);
   if(strlen(beta)>0)              printf("  BETA               <--  %s\n",beta);
   if(strlen(polynomial_degree)>0) printf("  POLYNOMIAL_DEGREE  <--  %s\n",polynomial_degree);
   if(strlen(shape_parameter)>0)   printf("  SHAPE_PARAMETER    <--  %s\n",shape_parameter);
   printf("  DATA_DIMENSION     <--  %s\n",data_dimension);   
   printf("  FILE_INPUT_CENTERS <--  %s\n",input_centers);   
   printf("  FILE_INPUT_DATA    <--  %s\n",input_data);
   printf("  FILE_OUTPUT_DATA   <--  %s\n",output_data);  
   printf("\n");       
   
      
//validate the token with the rbf selected
  if(strcmp(rbf,"MQ")==0)
      rbf_type = 1;
  else if(strcmp(rbf,"IMQ")==0)
      rbf_type = 2;
  else if(strcmp(rbf,"TPS")==0)
      rbf_type = 3;
  else if(strcmp(rbf,"POT")==0)
      rbf_type = 4;
  else if(strcmp(rbf,"GAU")==0)
      rbf_type = 5;
  else{
       fprintf(stderr,"ERROR: the rbf = %s, is not a valid name.\n\n",rbf);
       exit(1);
  }   
  
//validate the beta value 
  beta_factor = 1;  // to avoid compiler warning

  if(strlen(beta)>0){
     
      beta_factor = atoi(beta);
     
      if(beta_factor<0){
        fprintf(stderr,"ERROR: the value of beta = %d, must be a positive value.\n\n",beta_factor);
         exit(1);
      }    
     
      beta_band = 1;
  }   
  else{
     
       beta_band = 0; 
       
  }   
        
//validate the degree of the polynomial
    pol_degree = 0;  // to avoid compiler warning
  
  if(strlen(polynomial_degree)>0){
     
     pol_degree = atoi(polynomial_degree);
     
     if(pol_degree<0){
       fprintf(stderr,"ERROR: the polynomial_degree = %d, must be a positive value.\n\n",pol_degree);
       exit(1);
     }
     
     pol_degree_band = 1;
  
  }else{
     
      pol_degree_band  = 0;
      
  }      
   
//validate the value of the shape parameter
     c_value  = 1;  // to avoid compiler warning

   if(strlen(shape_parameter)>0){
     
     c_value = atof(shape_parameter);
    
     if(c_value<0.0){
          fprintf(stderr,"ERROR: the shape_parameter = %f, must be a positive value.\n\n",c_value);
          exit(1);
      }   
      
     c_value_band = 1; 
     
   }else{
     
     c_value_band = 0;
   
   }     

//validate the dimension of the data
   if(strlen(data_dimension)>0){
     
     dim = atoi(data_dimension);
    
     if(dim<1){
          fprintf(stderr,"ERROR: the data_dimension = %d, must be a positive value.\n\n",dim);
          exit(1);
      }   
     
   }else{
     
         fprintf(stderr,"ERROR: the data_dimension = '%s', not founded.\n\n",data_dimension);
          exit(1);    
   
   }     


//validate FILE_INPUT_CENTERS
   in = fopen(input_centers,"r");
   
   if(in == NULL){
     printf("\nERROR: unabled to open the filename '%s'.\n\n",input_centers);
     exit(1);   
   }
   
   fclose(in);
   
//validate FILE_INPUT_DATA
   in = fopen(input_data,"r");
   
   if(in == NULL){
     fprintf(stderr,"\nERROR: unabled to open the filename '%s'.\n\n",input_data);
     exit(1);   
   }
   
   fclose(in);   
   
//Now, we procede to interpolate the scattered data
  TPS<double> tps;
  POT<double> pot;   
  GAU<double> gau;
  IMQ<double> imq;
  MQ<double>   mq;
  
  
  switch(rbf_type)
  {
     case  1:
             interpolate( mq             , 
                        dim             ,
                        beta_band       , beta_factor,
                        pol_degree_band , pol_degree,
                        c_value_band    , c_value,
                        input_centers   ,
                        input_data      ,
                        output_data
                      );   
             break;        
      case  2:
             interpolate( imq , 
                        dim             ,
                        beta_band       , beta_factor,
                        pol_degree_band , pol_degree,
                        c_value_band    , c_value,
                        input_centers   ,
                        input_data      ,
                        output_data
                      );   
             break;      
      case  3:
             interpolate( tps , 
                        dim             ,
                        beta_band       , beta_factor,
                        pol_degree_band , pol_degree,
                        c_value_band    , c_value,
                        input_centers   ,
                        input_data      ,
                        output_data
                      );   
             break;      
      case  4:
             interpolate( pot , 
                        dim             ,
                        beta_band       , beta_factor,
                        pol_degree_band , pol_degree,
                        c_value_band    , c_value,
                        input_centers   ,
                        input_data      ,
                        output_data
                      );   
             break;      
      case  5:
             interpolate( gau , 
                        dim             ,
                        beta_band       , beta_factor,
                        pol_degree_band , pol_degree,
                        c_value_band    , c_value,
                        input_centers   ,
                        input_data      ,
                        output_data
                      );   
             break;      
  }


 return 0;
}

