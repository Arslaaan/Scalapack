#include <cstdio>
#include "mpi.h"
#include <complex>
#include <iostream>
#include "scalapack.h"
#include <iomanip>

using namespace std;

class DistributedMatrix{
  complex<double>* data;
  int rows,cols;//размер матрицы
  int procrows,proccols;//размер сетки процессов
  int myProcRow,myProcCol;
  int myProcRows,myProcCols;//количество строк и столбцов матрицы в процессе
  int myProcRowsOffset,myProcColsOffset;
  int myProcSize;
  int rank,nproc;
  int icon;
  int blocksize;
  int* desc;
public:
  DistributedMatrix(int rows,int cols,int lrank,int lnproc,int lblocksize,int licon) {
    rank = lrank;
    nproc = lnproc;
    blocksize = lblocksize;
    icon = licon;
    this->rows = rows;
    this->cols = cols;
    int divA;
    for (divA = (int) sqrt(nproc); (nproc % divA != 0); divA--);
    int divB = nproc / divA;
    if(rows > cols){
      procrows = divB;
      proccols = divA;
    }
    else {
      procrows = divA;
      proccols = divB;
    }
    Create();
  }

  void Create() {
    int x,y;
    int root = 0;
    Cblacs_gridinfo(icon, &x, &y, &myProcCol, &myProcRow);
    myProcRows = numroc_(&rows, &blocksize, &myProcRow, &root, &procrows);
    myProcCols = numroc_(&cols, &blocksize, &myProcCol, &root, &proccols);
    myProcRowsOffset = npreroc_(&rows, &blocksize, &myProcRow, &root, &procrows);
    myProcColsOffset = npreroc_(&cols, &blocksize, &myProcCol, &root, &proccols);//может в конце procrows
    myProcSize = myProcRows * myProcCols;
    //cout<<"myprocrows/cols "<<myProcRows<<" "<<myProcCols<<endl;
    data = new complex<double>[myProcSize];
    desc = descinit();
  }

  void Populate() {
    MPI_File thefile;
    MPI_Offset filesize;
    MPI_File_open(MPI_COMM_WORLD,"matrix",MPI_MODE_RDONLY,MPI_INFO_NULL,&thefile);
    MPI_File_get_size(thefile,&filesize);
  	for (int i = 0; i < myProcRows; i++) {
      double* tmp_buf = new double[2*myProcCols];
      MPI_File_set_view(thefile,2*sizeof(double)*(cols*(myProcRowsOffset+i)+myProcColsOffset),MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
      MPI_File_read(thefile,tmp_buf,2*myProcCols,MPI_DOUBLE,MPI_STATUS_IGNORE);
  		for (int j = 0; j < myProcCols; j++) {
  			data[i*myProcCols + j] = complex<double>(tmp_buf[j*2],0);
  		}
      delete(tmp_buf);
  	}
  }

  void Fill_diag() {
    if(myProcCol != myProcRow) {
      return;
    }
    for(int i = 0; i < myProcRows; i++) {
      data[i*(myProcCols+1)] = i+myProcRowsOffset+1;
    }
  }

  void Fill(complex<double> c) {
    for (int i = 0; i < myProcRows; i++) {
      for (int j = 0; j < myProcCols; j++) {
        data[i*myProcCols + j] = c;
      }
    }
  }

  void Print() {
    for(int i = 0; i < procrows; ++i) {
      for(int j = 0; j < proccols; ++j) {
        MPI_Barrier(MPI_COMM_WORLD);
        if(i*proccols + j == rank) {
          cout<<"coordinates of processes grid:"<<endl;
          cout<<"("<<myProcRow<<","<<myProcCol<<")"<<endl;
          for(int x = 0; x < myProcRows; ++x) {
            for(int y = 0; y < myProcCols; ++y) {
              cout<< "["<<data[x*myProcCols + y]<<"]";
            }
            cout<<endl;
          }
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  complex<double>* getData() {
    return data;
  }

  void setData(complex<double>* new_data) {
    for(int i = 0; i < myProcRows; i++) {
      for(int j = 0; j < myProcCols; j++) {
        int idx = i*myProcCols + j;
        data[idx] = new_data[idx + myProcRowsOffset*cols + myProcColsOffset];
      }
    }
  }

  int* getDesc() {
    return desc;
  }

  /*
   * result: C=C+AB
   * C изменится,
   * C меняется во время процедуры, так что не использовать ее в качестве A или B
   */
  void Multiply(DistributedMatrix& B, DistributedMatrix& C) {
    complex<double> cdOne(1.0,0.0);
    complex<double> cdZero(0.0, 0.0);
    int ione = 1;
    if(rank == 0)
      cout<<"--------------A--------------"<<endl;
    Print();
    if(rank == 0)
      cout<<"--------------B--------------"<<endl;
    B.Print();

    pzgemm_((char*)"N", (char*)"N", &rows, &rows, &rows,
            &cdOne, data, &ione, &ione, desc,
            B.getData(), &ione, &ione, B.getDesc(),
            &cdZero, C.getData(), &ione, &ione, C.getDesc());
  }

  int* descinit() {
    int* descme = new int[9];
    int iZero = 0;
    int info;
    int leadDim = max(myProcRows, myProcCols);
    //cout<<"myProcRows/Cols="<<myProcRows<<" "<<myProcCols<<endl;
    //cout<<"leadDim="<<leadDim<<endl;
    descinit_(descme, &rows, &cols, &blocksize, &blocksize, &iZero, &iZero, &icon, &leadDim, &info);
    if(info < 0)
      cout<<"DESC INIT ERROR\n";
    return descme;
  }

  void createSVD(DistributedMatrix& middle, DistributedMatrix& u, DistributedMatrix& vt) {
    // TODO test it

    double* s = new double [rows];//
    complex<double> singularValues[rows*cols];

    int iOne = 1;
    complex<double> *work = new complex<double>[1];
    int info = 0, lwork = -1;
    double *rwork = new double[1 + 4 * max(rows, cols)];

    pzgesvd_((char*) "V", (char*) "V",
      &rows, &cols, data, &iOne, &iOne, desc, s,
      u.getData(), &iOne, &iOne, u.getDesc(),
      vt.getData(), &iOne, &iOne, vt.getDesc(),
      work, &lwork, rwork, &info);

    lwork = (int) work[0].real();
    delete[] work;
    work = new complex<double>[lwork];

    pzgesvd_((char*) "V", (char*) "V",
      &rows, &cols, data, &iOne, &iOne, desc, s,
      u.getData(), &iOne, &iOne, u.getDesc(),
      vt.getData(), &iOne, &iOne, vt.getDesc(),
      work, &lwork, rwork, &info);

    //delete[] work, rwork;
  }


  ~DistributedMatrix() {
    //Cblacs_gridexit(icon);
    //Cblacs_exit(0);
  }

};

void test_mult(int rank, int nproc, int icon) {
  DistributedMatrix A(6,6,rank,nproc,3,icon);
  DistributedMatrix B(6,6,rank,nproc,3,icon);
  DistributedMatrix C(6,6,rank,nproc,3,icon);
  A.Populate();
  B.Populate();
  C.Fill(0);
  A.Multiply(B,C);
  C.Print();
}

void test_svd(int rank, int nproc, int icon) {
  DistributedMatrix A(6,6,rank,nproc,3,icon);
  DistributedMatrix B(6,6,rank,nproc,3,icon);
  DistributedMatrix C(6,6,rank,nproc,3,icon);
  DistributedMatrix S(6,6,rank,nproc,3,icon);
  A.Fill_diag();
  B.Fill(0);
  C.Fill(0);

  A.createSVD(S, B, C);

  A.Print();
  B.Print();
  C.Print();
  S.Print();
}


int main(int argc,char** argv){
  // init
  int rank,nproc;
  MPI_Init(&argc, &argv);
  Cblacs_pinfo(&rank, &nproc);
  int icon = -2;
  Cblacs_get(-1, 0, &icon);
  Cblacs_gridinit(&icon, (char *) "r", 2, 2);

  test_mult(rank, nproc, icon);
  //test_svd(rank, nproc, icon);

  // exit
  //Cblacs_gridexit(icon);
  MPI_Barrier(MPI_COMM_WORLD);
  Cblacs_exit(0);
}
