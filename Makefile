CC = g++ --std=c++11 -g
INC = -I/usr/include/eigen3

default:	lorenz95 lorenz63 libnilss.so

clean:
	rm -f *.o *.so lorenz63 lorenz95

# ------------------------ Lorenz63 demo program ----------------------------- #

lorenz63:	lorenz63.o libnilss.so
	gfortran -Wl,-R -Wl,. -o $@ $^

lorenz63.o:	lorenz63.f90
	gfortran -c -g $^

# ------------------------ Lorenz95 demo program ----------------------------- #

lorenz95:	lorenz95.o libnilss.so
	gfortran -Wl,-R -Wl,. -o $@ $^

lorenz95.o:	lorenz95.f90
	gfortran -c -g $^

# ------------------------ NILSS shared library ------------------------------ #

libnilss.so:	nilss_fortran.o nilss.o nilss_solver.o
	$(CC) -shared -o $@ $^

nilss_fortran.o:	nilss_fortran.cpp nilss.h
	$(CC) $(INC) -c -fPIC $<

nilss.o:	nilss.cpp nilss.h nilss_solver.h
	$(CC) $(INC) -c -fPIC $<

nilss_solver.o:	nilss_solver.cpp nilss_solver.h
	$(CC) $(INC) -c -fPIC $<
