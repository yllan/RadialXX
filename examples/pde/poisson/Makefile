include ../../Makefile.in

help: 
	@ echo ""
	@ echo "You have to type as follows"
	@ echo "          make <target>"
	@ echo "where <target> is one of:"
	@ echo " "
	@ echo "poisson   ---- Compile poisson problem"
	@ echo "poissonall  -- Compile the poisson_all.cpp problem"
	@ echo "all         -- Compile the previous examples"
	@ echo " "
	@ echo "clean   ------ Erase all the *.o, data and the executable files"
	@ echo "cleanobj ----- Erase all the *.o"
	@ echo " "

all: poisson poissonall  

poisson: poisson.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -c poisson.cpp
	$(CXX)   poisson.o -o poisson

poissonall: poisson_all.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -c poisson_all.cpp
	$(CXX)   poisson_all.o -o poissonall

clean: cleanobj; rm -f poisson  poissonall data
 
cleanobj:
	rm -f *.o 
