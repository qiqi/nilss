CTG = cxxtestgen --error-printer
CXX = g++ --std=c++11 -g
INC = -I/usr/include/eigen3 -I..

default:	unittests testFortran

clean:
	rm -f testLSS.c *.o unittest

unittests:	testLSS.c ../libnilss.so
	$(CXX) $(INC) -Wl,-R -Wl,.. -o unittest testLSS.c -L.. -lnilss

testLSS.c:	testLSS.h
	$(CTG) -o testLSS.c testLSS.h

testFortran:	testFortran.f90 testFortran_c.o
	gfortran -Wl,-R -Wl,. -o testFortran $^ -lstdc++

testFortran_c.o:	testFortran_c.cpp
	$(CXX) -c $<
