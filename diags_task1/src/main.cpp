#include <cstdio>
#include "mpi.h"
#include <complex>
#include <iostream>
#include "scalapack.h"
#include <iomanip>

using namespace std;

class Matrix{
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
  Matrix(int rows,int cols,int lrank,int lnproc,int lblocksize,int licon){
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
    else{
      procrows = divA;
      proccols = divB;
    }
    //cout<<"rows="<<procrows<<" cols="<<proccols<<endl;
    Create();
    /*
    for(int i =0;i<9;++i)
      if(rank == 1)
        cout<<desc[i]<<" ";
    */
  }

  void Create(){
    int x,y;
    int root = 0;
    Cblacs_gridinit(&icon, (char *) "Column", procrows, proccols);
    Cblacs_gridinfo(icon, &x, &y, &myProcRow, &myProcCol);
    myProcRows = numroc_(&rows, &blocksize, &myProcRow, &root, &procrows);
    myProcCols = numroc_(&cols, &blocksize, &myProcCol, &root, &proccols);
    myProcRowsOffset = npreroc_(&rows, &blocksize, &myProcRow, &root, &procrows);
    myProcColsOffset = npreroc_(&cols, &blocksize, &myProcCol, &root, &proccols);//может в конце procrows
    myProcSize = myProcRows * myProcCols;
    //cout<<"myprocrows/cols "<<myProcRows<<" "<<myProcCols<<endl;
    data = new complex<double>[myProcSize];
    desc = descinit();
  }

  void Populate(){
    MPI_File thefile;
    MPI_Offset filesize;
    MPI_File_open(MPI_COMM_WORLD,"matrix",MPI_MODE_RDONLY,MPI_INFO_NULL,&thefile);
    MPI_File_get_size(thefile,&filesize);
  	for (int i = 0; i < myProcRows; i++) {
      double* tmp_buf = new double[2*myProcCols];
      MPI_File_set_view(thefile,2*sizeof(double)*(cols*(myProcRowsOffset+i)+myProcColsOffset),MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
      MPI_File_read(thefile,tmp_buf,2*myProcCols,MPI_DOUBLE,MPI_STATUS_IGNORE);
  		for (int j = 0; j < myProcCols; j++) {
  			data[i*myProcCols + j] = complex<double>(tmp_buf[j*2],tmp_buf[j*2 + 1]);
  		}
      delete(tmp_buf);
  	}
  }

  void Print(){
    for(int i = 0;i < procrows;++i)
      for(int j = 0;j < proccols;++j){
        MPI_Barrier(MPI_COMM_WORLD);
        if(i*proccols + j == rank){
          cout<<"coordinates of processes grid:"<<endl;
          cout<<"("<<myProcRow<<","<<myProcCol<<")"<<endl;
          for(int x = 0;x < myProcRows;++x){
            for(int y = 0;y < myProcCols;++y){
              cout<< "["<<data[x*myProcCols + y]<<"]";
            }
            cout<<endl;
          }
        }
      }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  complex<double>* getData(){
    return data;
  }

  int* getDesc(){
    return desc;
  }

  void Multiply(Matrix& B, Matrix& C){
    complex<double> cdOne(1.0,0.0);
    int ione = 1;

    pzgemm_((char*)"N", (char*)"N", &rows, &rows, &rows, &cdOne, data, &ione, &ione, desc, B.getData(), &ione, &ione, B.getDesc(), &cdOne, C.getData(), &ione, &ione, C.getDesc());
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
/*
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
*/

  ~Matrix(){
    //Cblacs_gridexit(icon);
    //Cblacs_exit(0);
  }

};


int main(int argc,char** argv){
  int rank,nproc;
  MPI_Init(&argc, &argv);
  Cblacs_pinfo(&rank, &nproc);
  int icon = -2;
  Cblacs_get(-1, 0, &icon);
  Matrix A(6,6,rank,nproc,3,icon);
  Matrix B(6,6,rank,nproc,3,icon);
  Matrix C(6,6,rank,nproc,3,icon);
  A.Populate();
  B.Populate();
  //A.Print();
  //B.Print();
  A.Multiply(B,C);
  //C.Print();
  //Cblacs_gridexit(icon);
  Cblacs_exit(0);
}
