CTG = cxxtestgen --error-printer
CXX = g++ --std=c++11 -g
INC = -I/usr/include/eigen3 -I..

default:	unittests testFortran

clean:
	rm -f test_*.c *.o unittest

unittests:	test_nilss.c ../libnilss.so
	$(CXX) $(INC) -Wl,-R -Wl,.. -o unittests test_nilss.c -L.. -lnilss

test_nilss.c:	test_nilss.h
	$(CTG) -o test_nilss.c test_nilss.h

testFortran:	testFortran.f90 testFortran_c.o
	gfortran -Wl,-R -Wl,. -o testFortran $^ -lstdc++

testFortran_c.o:	testFortran_c.cpp
	$(CXX) -c $<
