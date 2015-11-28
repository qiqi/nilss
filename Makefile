CC = g++ --std=c++11 -g
INC = -I/usr/include/eigen3

default:	lib lib/libnilss.so

clean:
	rm -f *.o lib/*.so

lib:
	mkdir lib

lib/libnilss.so:	nilss_fortran.o nilss.o nilss_solver.o
	$(CC) -shared -o $@ $^

nilss_fortran.o:	nilss_fortran.cpp nilss.h
	$(CC) $(INC) -c -fPIC $<

nilss.o:	nilss.cpp nilss.h nilss_solver.h
	$(CC) $(INC) -c -fPIC $<

nilss_solver.o:	nilss_solver.cpp nilss_solver.h
	$(CC) $(INC) -c -fPIC $<
