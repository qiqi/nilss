CTG = cxxtestgen --error-printer
CXX = g++ --std=c++11 -g
INC = -I/usr/include/eigen3 -I..

default:	unittests testFortran

clean:
	rm -f _unittests*.c *.o unittests testFortran

unittests:	_unittests.c ../lib/libnilss.so
	$(CXX) $(INC) -Wl,-R -Wl,../lib -o $@ $^

_unittests.c:	test_nilss.h test_nilss_solver.h
	$(CTG) -o $@ $^

testFortran:	testFortran.f90 testFortran_c.o
	gfortran -Wl,-R -Wl,. -o $@ $^ -lstdc++

testFortran_c.o:	testFortran_c.cpp
	$(CXX) -c $<
