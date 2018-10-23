OPT=	-O2 
FLAGS=	-DUSE_COMPLEX -DUSE_DOUBLE -g

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


all: bin

bin: bin/main

#lib: lib/smatrix.a lib/soutput.a

build/main.o: src/main.cpp Makefile
	mkdir -p build
	mpicxx -o $@ $(OPT) -c src/main.cpp $(FLAGS)

bin/main: build/main.o Makefile
	mkdir -p bin
	mpicxx -o $@ $(OPT) build/main.o $(LIBS)

test: bin/main
	mpiexec -np 4 ./bin/main

clean:
	rm -rf build
	find . -name \*~ -delete

clear:
	rm -rf bin lib doc
