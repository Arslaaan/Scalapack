OPT=	-O2
FLAGS=	-DUSE_COMPLEX -DUSE_DOUBLE -g
N=6
dt=0.5
iterations=5
nproc=4

ifneq (,$(findstring lomonosov,$(shell hostname)))
	LIBS=	-lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_blacs_intelmpi_lp64 \
		-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
	EFLAGS=
else
	LIBPATH=/usr/local/lib/
	LIBS=	$(LIBPATH)libscalapack.a \
		-lpthread -lblas -llapack -lm -lgfortran -fPIC
	EFLAGS=
endif

matrix:
	g++ -Wall create_matrix.cpp -o create_matrix
	./create_matrix $(N) $(N) matrix line
	./create_matrix $(N) $(N) hamilton max_coord
	./create_matrix $(N) $(N) rho vector

all: bin

bin: lib bin/main

lib:
	mkdir -p lib
	ar rcs lib/DistributedMatrix.a build/main.o

build/main.o: src/main.cpp Makefile
	mkdir -p build
	mpicxx -o $@ $(OPT) -c src/main.cpp $(FLAGS)

bin/main: build/main.o Makefile
	mkdir -p bin
	mpicxx -o $@ $(OPT) build/main.o lib/DistributedMatrix.a $(LIBS)

test: lib bin/main src/main.cpp
	mpiexec -np $(nproc) ./bin/main $(dt) $(iterations) $(N)

memcheck: bin/main
	mpiexec -np $(nproc) valgrind ./bin/main > log.txt

clean:
	rm -rf build
	find . -name \*~ -delete

clear:
	rm -rf bin lib doc
	rm -f matrix create_matrix hamilton rho
