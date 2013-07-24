include ../../Makefile.in

help: 
	@ echo ""
	@ echo "You have to type as follows"
	@ echo "          make <target>"
	@ echo "where <target> is one of:"
	@ echo " "
	@ echo "burger    ---- Compile burger problem"
	@ echo "bur_blas  ---- Compile burger with BLAS and LAPACK"
	@ echo "all       ---- Compile the previous examples"
	@ echo " "
	@ echo "clean   ------ Erase all the *.o, data* and the executable files"
	@ echo "cleanobj ----- Erase all the *.o"
	@ echo " "

all:  burger  bur_blas

burger: burger.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -c burger.cpp
	$(CXX)   burger.o -o burger
	
bur_blas: burger.cpp
	$(CXX)  $(CXXFLAGS)  -I$(INCLUDES) -DWITH_LAPACK -c burger.cpp
	$(CXX) -lblas -llapack burger.o  -o bur_blas 
	
clean: cleanobj; rm -f burger; rm -f bur_blas; rm  -f data/data*
 
cleanobj:
	rm -f *.o 
