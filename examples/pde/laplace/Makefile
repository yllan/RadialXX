include ../../Makefile.in

help: 
	@ echo ""
	@ echo "You have to type as follows"
	@ echo "          make <target>"
	@ echo "where <target> is one of:"
	@ echo " "
	@ echo "laplace   ---- Compile poisson problem"
	@ echo " "
	@ echo "clean   ------ Erase all the *.o, data and the executable files"
	@ echo "cleanobj ----- Erase all the *.o"
	@ echo " "

all: laplace

laplace: laplace.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -c laplace.cpp
	$(CXX)   laplace.o -o laplace

clean: cleanobj; rm -f laplace; rm -f data
 
cleanobj:
	rm -f *.o 
