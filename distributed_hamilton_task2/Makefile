OPT=	-O2
FLAGS=	-DUSE_COMPLEX -DUSE_DOUBLE -g
nproc = 4

ifneq (,$(findstring lomonosov,$(shell hostname)))
	LIBS=	-lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_blacs_intelmpi_lp64 \
		-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
	EFLAGS=
else
	LIBPATH=/usr/lib/
	LIBS=	$(LIBPATH)libscalapack.so \
		-lpthread -lblas -llapack -lm -lgfortran -fPIC
	EFLAGS=
endif

all: bin/main

build/main.o: src/main.cpp
	mkdir -p build
	mpicxx -o $@ $(OPT) -c src/main.cpp $(FLAGS)

build/DistributedMatrix.o:src/DistributedMatrix.cpp src/DistributedMatrix.h
	mpicxx -o $@ $(OPT) -c src/DistributedMatrix.cpp $(FLAGS)

build/KetVector.o:src/KetVector.cpp src/DistributedMatrix.h
	mpicxx -o $@ $(OPT) -c src/KetVector.cpp $(FLAGS)

bin/main: build/main.o build/DistributedMatrix.o build/KetVector.o
	mkdir -p bin
	mpicxx -o $@ $(OPT) $^ $(LIBS)

test: bin/main src/main.cpp
	mpiexec -np $(nproc) --use-hwthread-cpus ./bin/main

memcheck: bin/main
	mpiexec -np $(nproc) valgrind ./bin/main > log.txt

clear:
	rm -rf bin build
	find . -name \*~ -delete
