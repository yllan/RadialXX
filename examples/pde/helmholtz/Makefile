include ../../Makefile.in

help: 
	@ echo ""
	@ echo "You have to type as follows"
	@ echo "          make <target>"
	@ echo "where <target> is one of:"
	@ echo " "
	@ echo "helm      ---- Compile helmholtz problem"
	@ echo " "
	@ echo "clean   ------ Erase all the *.o, data and the executable files"
	@ echo "cleanobj ----- Erase all the *.o"
	@ echo " "

helm: helmholtz.cpp
	$(CXX)  $(CXXFLAGS) -I$(INCLUDES) -c helmholtz.cpp
	$(CXX)   helmholtz.o -o helmholtz

clean: cleanobj; rm -f helmholtz  data
 
cleanobj:
	rm -f *.o 
