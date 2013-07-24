include ../../Makefile.in

help: 
	@ echo ""
	@ echo "You have to type as follows"
	@ echo "          make <target>"
	@ echo "where <target> is one of:"
	@ echo " "
	@ echo "eje1    ----- Compile eje-2d-mq problem"
	@ echo "eje1g   ----- Compile eje-2d-mq problem wit gnuplot output"
	@ echo "eje2    ----- Compile eje-2d-tps problem"	
	@ echo "ejeall  ----- Compile eje-2d-all problem"
	@ echo "all     ----- Compile the previous examples"
	@ echo " "
	@ echo "clean   ----- Erase all the *.o, output data and the executable files"
	@ echo "cleanobj ---- Erase all the *.o"
	@ echo " "

all: eje1 eje1g eje2 ejeall

eje1: eje-2d-mq.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -c eje-2d-mq.cpp
	$(CXX)   eje-2d-mq.o -o eje1
	
eje1g: eje-2d-mq.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -DWITH_GNUPLOT -c eje-2d-mq.cpp
	$(CXX)   eje-2d-mq.o -o eje1g
	
eje2: eje-2d-tps.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -c eje-2d-tps.cpp
	$(CXX)   eje-2d-tps.o -o eje2	
	
ejeall: eje-2d-all.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -c eje-2d-all.cpp
	$(CXX)   eje-2d-all.o -o ejeall	

clean: cleanobj; rm -f eje1 eje2  eje1g ejeall
 
cleanobj:
	rm -f *.o; rm -f *.dat 
