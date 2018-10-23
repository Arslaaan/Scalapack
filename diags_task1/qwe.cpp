#include <cstdio>

class Matrix {
private:
  int rows, cols; // size of matrix
  int procRows, procCols; // size of processes grid
  int myProcRow, myProcCol;
  int rank, nproc;
  int icon;
  int blockSize;
  int *desc;
  complex<double> *data;
  complex<double> *dataCopy;
public:
  descinit(complex<double>);
  Matrix(int rows, int cols, int lrank, int lnproc, int lblockSize) {
    rank = lrank;
    nproc = lnproc;
    blockSize = lblockSize;
    data = new complex<double>[rows*cols];
    desc = descinit(data);

    int divA;
    for (divA = (int) sqrt(nproc); (nproc % divA != 0); divA--) {
      // empty cycle
    }
    int divB = nproc / divA;
    if(rows > cols) {
      procRows = divB;
      procCols = divA;
    }
    else {
      procRows = divA;
      procCols = divB;
    }
    Cblacs_get(-1, 0, &icon);
    Cblacs_gridinit(&icon, (char *) "Column", procRows, procCols);
  }

  void set() {
    dataCopy = new complex<double>[procDim];
    for(int col = 0; col < procCols; col++) {
      for(int row = 0; row < procRows; row++) {
        // TODO test it
        dataCopy[procRows * col + row] = data[(col+procColsOffset)*rows + (row+procRowsOffset)*(cols-procCols)];
      }
    }
  }

  int *descinit(complex<double> me) {
    int descme[9];
    int iZero = 0;
    int info;
    int leadDim = max(rows, cols);
    descinit_(descme, &rows, &cols, &blockSize, &blockSize, &iZero, &iZero, &rank, &leadDim, &info);
    return descme;
  }

  void createSVD(double *s, complex<double> *u, complex<double> *vt) {
    // TODO test it

    int *descu = descinit(u), *descvt = descinit(vt);

    int iZero = 0;
    complex<double> *work = new complex<double>[1];
    int info = 0, lwork = -1;
    double *rwork = new double[1 + 4 * max(rows, cols)];

    pzgesvd_(char* "V", char* "V", 
      &rows, &cols, dataCopy, &iZero, &iZero, desc, s, 
      u, &iZero, &iZero, descu,
      vt, &iZero, &iZero, descvt,
      work, lwork, rwork, &info);
  }

}

int main(int argc,char** argv) {
  int rank,nproc;
  MPI_Init(&argc, &argv);
  Cblacs_pinfo(&rank, &nproc);

}
