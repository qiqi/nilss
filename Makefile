CC = g++ --std=c++11

default:	lorenz63 libnilss.so

clean:
	rm -f *.o *.so lorenz63

# ------------------------ Lorenz63 demo program ----------------------------- #

lorenz63:	libnilss.so lorenz63.o
	gfortran -Wl,-R -Wl,. -o lorenz63 lorenz63.o -L. -lnilss

lorenz63.o:	lorenz63.f90
	gfortran -c lorenz63.f90

# ------------------------ NILSS shared library ------------------------------ #

libnilss.so:	nilss_fortran.o nilss.o
	$(CC) -shared -o libnilss.so nilss.o nilss_fortran.o

nilss_fortran.o:	nilss_fortran.cpp nilss.h
	$(CC) -c -fPIC nilss_fortran.cpp

nilss.o:	nilss.cpp nilss.h
	$(CC) -I/usr/include/eigen3 -c -fPIC nilss.cpp
