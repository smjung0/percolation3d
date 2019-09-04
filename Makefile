CC = icpc
OPT_FLAGS =-O2 -xHost -qopenmp

perc3d : perc3d_main.o perc3d_solve.o perc3d_structure.o perc3d.o cnt.o
	icpc -std=c++11 ${OPT_FLAGS} -o perc3d perc3d_main.o perc3d_solve.o perc3d_structure.o perc3d.o cnt.o -liomp5 -lpthread -lm -ldl

perc3d_main.o : perc3d_main.cpp
	icpc -std=c++11 ${OPT_FLAGS} -c -o perc3d_main.o perc3d_main.cpp -liomp5 -lpthread -lm -ldl 

perc3d_solve.o : perc3d_solve.cpp
	icpc -std=c++11 ${OPT_FLAGS} -c -o perc3d_solve.o perc3d_solve.cpp -liomp5 -lpthread -lm -ldl

perc3d_structure.o : perc3d_structure.cpp
	icpc -std=c++11 ${OPT_FLAGS} -c -o perc3d_structure.o perc3d_structure.cpp -liomp5 -lpthread -lm -ldl

perc3d.o : perc3d.cpp
	icpc -std=c++11 ${OPT_FLAGS} -c -o perc3d.o perc3d.cpp -liomp5 -lpthread -lm -ldl 

cnt.o : cnt.cpp 
	icpc -std=c++11 ${OPT_FLAGS} -c -o cnt.o cnt.cpp -liomp5 -lpthread -lm -ldl 

clean :
	rm *.o



