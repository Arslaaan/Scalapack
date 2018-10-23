N=6
matrix:
	rm -f matrix
	g++ -Wall create_matrix.cpp -o create_matrix
	./create_matrix $N $N matrix
compile:
	mpicxx -o main -Wall -c  -O2 main.cpp -DUSE_COMPLEX -DUSE_DOUBLE -L/usr/local/lib/scalapack.a  
	mpiexec -np 4 ./main
