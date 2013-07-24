include ../../Makefile.in

help: 
	@ echo ""
	@ echo "You have to type as follows"
	@ echo "          make <target>"
	@ echo "where <target> is one of:"
	@ echo " "
	@ echo "eje1    ----- Compile eje-1 problem"
	@ echo " "
	@ echo "clean   ----- Erase all the *.o, data files and the executable files"
	@ echo "cleanobj ---- Erase all the *.o"
	@ echo " "

eje1: eje-1.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -c eje-1.cpp
	$(CXX)  eje-1.o   -o eje1 
	

cleanobj:
	rm -f *.o 

clean:
	rm -f *.o *.dat eje1
