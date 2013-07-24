include ../../Makefile.in

help: 
	@ echo ""
	@ echo "You have to type as follows"
	@ echo "          make <target>"
	@ echo "where <target> is one of:"
	@ echo " "
	@ echo "circle    ---- Compile circle problem"
	@ echo " "
	@ echo "clean   ------ Erase all the *.o, data and the executable files"
	@ echo "cleanobj ----- Erase all the *.o"
	@ echo " "

all: circle

circle: circle.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -c circle.cpp
	$(CXX)   circle.o -o circle	

clean: cleanobj; rm -f circle  data
 
cleanobj:
	rm -f *.o 
