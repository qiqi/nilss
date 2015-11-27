default:	lorenz95 lorenz63

clean:
	rm -f *.o *.so lorenz63 lorenz95

lorenz95:	lorenz95.f90
	gfortran -Wl,-R -Wl,. -o $@ $^

lorenz63:	lorenz63.f90
	gfortran -Wl,-R -Wl,. -o $@ $^
