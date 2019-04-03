# Compilers
cc=gcc-8
cxx=g++-8

# Compilation flags
cflags=-Wall -Wextra -O3 -std=c++17 -fopenmp -march=native -msse -msse2 -mavx
lflags=-lfftw3 -lgsl -lgslcblas -lm

# Lists of files to be built
objs=
src=$(patsubst %.o,%.cpp,$(objs))
execs=test_nufftw fftw_benchmark

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm $(objs) $(execs)

%.o: %.cpp
	$(cxx) $(cflags) -c $<

test_nufftw: test_nufftw.cpp $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lflags)

fftw_benchmark: fftw_benchmark.cpp $(objs)
	$(cxx) $(cflags) -o $@ $^ $(lflags)

.PHONY: clean depend
