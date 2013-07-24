include ../../Makefile.in

help: 
	@ echo " "
	@ echo "You have to type as follows"
	@ echo "          make <target>"
	@ echo "where <target> is one of:"
	@ echo " "
	@ echo "advection ---- Compile advection problem"
	@ echo "adv_blas  ---- Compile advection with BLAS and LAPACK"
	@ echo "adv_inv    --- Compile advection using inv and BLAS and LAPACK"
	@ echo "all      ----- Compile previous examples"
	@ echo " "
	@ echo "clean   ------ Erase all the *.o, data* and the executable files"
	@ echo "cleanobj ----- Erase all the *.o"
	@ echo " "

advection: advection.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -c advection.cpp
	$(CXX)   advection.o -o advection
	
adv_blas: advection.cpp
	$(CXX)  $(CXXFLAGS)  -I$(INCLUDES) -DWITH_LAPACK -c advection.cpp
	$(CXX) -lblas -llapack advection.o -o advection_blas 
		
adv_inv: advection_inv.cpp
	$(CXX)  $(CXXFLAGS)  -I$(INCLUDES) -DWITH_LAPACK -c advection_inv.cpp
	$(CXX) -lblas -llapack advection_inv.o  -o advection_inv
	
all: advection  adv_blas  adv_inv

clean: cleanobj; rm -f advection; rm -f advection_blas; rm -f advection_inv; rm  -f data/data*
 
cleanobj:
	rm -f *.o 
	
