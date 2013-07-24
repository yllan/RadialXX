include ../../Makefile.in

help: 
	@ echo ""
	@ echo "You have to type as follows"
	@ echo "          make <target>"
	@ echo "where <target> is one of:"
	@ echo " "
	@ echo "interpola     ----- Compile interpolate problem"
	@ echo "interpola_blas  --- Compile interpolate with BLAS and LAPACK"
	@ echo " "
	@ echo "clean         ----- Erase all the *.o and the executable files"
	@ echo "cleanobj      ----- Erase all the *.o"
	@ echo " "

interpola: interpolate.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -c interpolate.cpp
	$(CXX)   interpolate.o -o interpolate

interpola_blas: interpolate.cpp
	$(CXX)  $(CXXFLAGS)  -I$(INCLUDES) -DWITH_LAPACK -c interpolate.cpp
	$(CXX) -lblas -llapack interpolate.o  -o interpolate_blas 
	
cleanobj:
	rm -f *.o 

clean:
	rm -f *.o 
	rm -f interpolate interpolate_blas
