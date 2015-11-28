CC = g++ --std=c++11 -g
INC = -I/usr/include/eigen3

default:	lorenz95 lorenz63

clean:
	rm -f *.o lorenz63 lorenz95

# ------------------------ Lorenz63 demo program ----------------------------- #

lorenz63:	lorenz63.o ../lib/libnilss.so
	gfortran -Wl,-R -Wl,. -o $@ $^

lorenz63.o:	lorenz63.f90
	gfortran -c -g $^

# ------------------------ Lorenz95 demo program ----------------------------- #

lorenz95:	lorenz95.o ../lib/libnilss.so
	gfortran -Wl,-R -Wl,. -o $@ $^

lorenz95.o:	lorenz95.f90
	gfortran -c -g $^

