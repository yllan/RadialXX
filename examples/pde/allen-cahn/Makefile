include ../../Makefile.in

help: 
	@ echo ""
	@ echo "You have to type as follows"
	@ echo "          make <target>"
	@ echo "where <target> is one of:"
	@ echo " "
	@ echo "allen     ---- Compile allen-cahn problem"
	@ echo "allen_blas  -- Compile allen-cahn with BLAS and LAPACK"
	@ echo " "
	@ echo "clean   ------ Erase all the *.o, data* and the executable files"
	@ echo "cleanobj ----- Erase all the *.o"
	@ echo " "

allen: allen-cahn.cpp
	$(CXX)  $(CXXFLAGS)  -I$(INCLUDES) -c allen-cahn.cpp
	$(CXX)  allen-cahn.o  -o allen
	
allen_blas: allen-cahn.cpp
	$(CXX)  $(CXXFLAGS)  -I$(INCLUDES) -DWITH_LAPACK -c allen-cahn.cpp
	$(CXX) -lblas -llapack allen-cahn.o  -o allen_blas 

clean: cleanobj; rm -f allen; rm -f allen_blas; rm  -f data/data*
 
cleanobj:
	rm -f *.o 
