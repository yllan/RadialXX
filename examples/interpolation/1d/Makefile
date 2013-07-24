include ../../Makefile.in

help: 
	@ echo ""
	@ echo "You have to type as follows"
	@ echo "          make <target>"
	@ echo "where <target> is one of:"
	@ echo " "
	@ echo "eje1    ----- Compile eje-1 problem"
	@ echo "eje1g   ----- Compiler eje-1 problem with output in gnuplot"
	@ echo "eje2     ---- Compile  eje-2 problem"
	@ echo "eje2g   ----- Compiler eje-2 problem with output in gnuplot"
	@ echo "all     ----- Compile the previous examples"
	@ echo " "
	@ echo "clean   ----- Erase all the *.o, output data and the executable files"
	@ echo "cleanobj ---- Erase all the *.o"
	@ echo " "

all: eje1 eje2 eje1g eje2g

eje1: eje-1.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -c eje-1.cpp
	$(CXX)   eje-1.o -o eje1
	
eje1g: eje-1.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -DWITH_GNUPLOT -c eje-1.cpp
	$(CXX)   eje-1.o -o eje1g	
	
eje2: eje-2.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -c eje-2.cpp
	$(CXX)   eje-2.o -o eje2	

eje2g: eje-2.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -DWITH_GNUPLOT -c eje-2.cpp
	$(CXX)   eje-2.o -o eje2g	

clean: cleanobj; rm -f eje1 eje2 eje1g eje2g
 
cleanobj:
	rm -f *.o; rm -f *.dat 
