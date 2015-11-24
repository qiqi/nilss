CC = g++ --std=c++11 -g
INC = -I/usr/include/eigen3

default:	lorenz63 libnilss.so

clean:
	rm -f *.o *.so lorenz63

# ------------------------ Lorenz63 demo program ----------------------------- #

lorenz63:	libnilss.so lorenz63.o
	gfortran -Wl,-R -Wl,. -o lorenz63 lorenz63.o -L. -lnilss

lorenz63.o:	lorenz63.f90
	gfortran -c lorenz63.f90

# ------------------------ NILSS shared library ------------------------------ #

libnilss.so:	nilss_fortran.o nilss.o nilss_solver.o
	$(CC) -shared -o libnilss.so nilss_fortran.o nilss.o nilss_solver.o

nilss_fortran.o:	nilss_fortran.cpp nilss.h
	$(CC) $(INC) -c -fPIC nilss_fortran.cpp

nilss.o:	nilss.cpp nilss.h nilss_solver.h
	$(CC) $(INC) -c -fPIC nilss.cpp

nilss_solver.o:	nilss_solver.cpp nilss_solver.h
	$(CC) $(INC) -c -fPIC nilss_solver.cpp
