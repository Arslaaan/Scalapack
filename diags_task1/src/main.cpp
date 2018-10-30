#include <cstdio>
#include "mpi.h"
#include <complex>
#include <iostream>
#include "scalapack.h"
#include <iomanip>
#include <string>

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

  void Fill(string name) {
    MPI_File thefile;
    MPI_Offset filesize;
    MPI_File_open(MPI_COMM_WORLD,name.c_str(),MPI_MODE_RDONLY,MPI_INFO_NULL,&thefile);
    MPI_File_get_size(thefile,&filesize);
    filesize = filesize/sizeof(double);
    const int matrix = rows*cols*2,qvector = rows*2,number = 2;
    if(filesize == matrix){
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
    else{
      if(filesize == qvector){
      }
      else{
        cout<<filesize<<endl;
        cout<<"err1"<<endl;
      }
    }
  }

  void Fill(complex<double> c) {
    for (int i = 0; i < myProcRows; i++) {
      for (int j = 0; j < myProcCols; j++) {
        data[i*myProcCols + j] = c;
      }
    }
  }

  /*
   * заполняет диагональ матрицы значениями из double*
   */
  void Fill(double* x){
    for (int i = 0; i < myProcRows; i++) {
  		for (int j = 0; j < myProcCols; j++) {
        if(myProcRow == myProcCol) {
          if(i == j)
    			   data[i*myProcCols + j] = x[i + myProcRow*blocksize];
        }
        else{
          data[i*myProcCols + j] = 0;
        }
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
          MPI_Barrier(MPI_COMM_WORLD);
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
      cout<<"--------------------------------"<<endl;
  }

  void Print(string key) {
    MPI_Barrier(MPI_COMM_WORLD);
    if(key == "diagonal"){
      for(int i = 0; i < nproc; ++i) {
        MPI_Barrier(MPI_COMM_WORLD);
        if((myProcRow == myProcCol) && (i == myProcRow)){
          for(int x = 0; x < myProcRows; ++x)
            if((myProcRow == procrows - 1) && (x == myProcRows - 1))
              cout<<data[x*myProcCols + x]<<flush;
            else
              cout<<data[x*myProcCols + x]<<","<<flush;
        }
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == nproc - 1)
      cout<<endl;
  }

  complex<double>* getData() {
    return data;
  }

  int* getDesc() {
    return desc;
  }

  /*
   * result: C=C+this*B
   * C изменится,
   * C меняется во время процедуры, так что не использовать ее в качестве this или B
   */
  void Multiply(DistributedMatrix& B, DistributedMatrix& C) {
    complex<double> cdOne(1.0,0.0);
    complex<double> cdZero(0.0, 0.0);
    int ione = 1;

    pzgemm_((char*)"N", (char*)"N", &rows, &rows, &rows,
            &cdOne, data, &ione, &ione, desc,
            B.getData(), &ione, &ione, B.getDesc(),
            &cdZero, C.getData(), &ione, &ione, C.getDesc());
  }

  void Multiply(DistributedMatrix& A,DistributedMatrix& B,DistributedMatrix& C){
    DistributedMatrix tmp1(6,6,rank,nproc,3,icon);
    DistributedMatrix tmp2(6,6,rank,nproc,3,icon);
    B.Multiply(C,tmp1);
    A.Multiply(tmp1,tmp2);
    *this = tmp2;
  }

  void conjTrans(DistributedMatrix& B) {
    B.Fill(static_cast<double>(0));
    complex<double> cdOne(1.0,0.0);
    complex<double> cdZero(0.0, 0.0);
    int ione = 1;
    DistributedMatrix I(rows,cols,rank,nproc,3,icon);
    double* ones = new double[rows];
    for(int i = 0;i < rows;++i)
      ones[i] = 1;
    I.Fill(ones);
    pzgemm_((char*)"C", (char*)"N", &rows, &rows, &rows,
            &cdOne, data, &ione, &ione, desc,
            I.getData(), &ione, &ione, I.getDesc(),
            &cdZero, B.getData(), &ione, &ione, B.getDesc());
/*
*/
  }


  void CreateSVD(DistributedMatrix& middle, DistributedMatrix& u, DistributedMatrix& vt) {
    // TODO test it

    double* s = new double [rows];

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

    middle.Fill(s);
    delete(s);
    //delete[] work, rwork;
  }

  void Exp(double dt){
    for(int i = 0;i < myProcRows;++i){
      data[i*myProcCols + i] = complex<double>(0,-1)*exp(data[i*myProcCols + i])*dt;
    }
  }

  DistributedMatrix operator = (DistributedMatrix& t){
    for(int i = 0;i < myProcRows;++i)
      for(int j = 0;j < myProcCols;++j){
        data[i*myProcCols + j] = t.getData()[i*myProcCols + j];
      }
  }

};

void test_svd(int rank, int nproc, int icon) {
  double dt = 0.5;
  int n = 4;
  DistributedMatrix H(6,6,rank,nproc,3,icon);
  DistributedMatrix rho(6,6,rank,nproc,3,icon);
  DistributedMatrix tmp(6,6,rank,nproc,3,icon);
  DistributedMatrix left(6,6,rank,nproc,3,icon);
  DistributedMatrix right(6,6,rank,nproc,3,icon);
  DistributedMatrix middle(6,6,rank,nproc,3,icon);
  DistributedMatrix unitary(6,6,rank,nproc,3,icon);
  DistributedMatrix transposed(6,6,rank,nproc,3,icon);
  H.Fill("hamilton");
  rho.Fill("matrix");
  left.Fill(static_cast<double>(0));
  middle.Fill(static_cast<double>(0));
  right.Fill(static_cast<double>(0));
  tmp.Fill(static_cast<double>(0));

  H.CreateSVD(middle, left, right);// H == left*middle*right
  middle.Exp(dt);
  unitary.Multiply(left,middle,right);//unitary = left*middle*right;

  //все выше верно работает
  unitary.conjTrans(transposed);
  transposed.Multiply(unitary,tmp);
  //должно печатать единичную матрицу,но не печатает,
  //ошибка где то в транспонировании
  tmp.Print();

/*
  for(int i = 0;i < n;++i){
    //rho.Print();
    rho.Print("diagonal");
    tmp.Multiply(transposed,rho,unitary);
    rho = tmp;
  }
*/
}


int main(int argc,char** argv){
  // init
  int rank,nproc;
  int gridrow = 2,gridcol = 2;
  MPI_Init(&argc, &argv);
  Cblacs_pinfo(&rank, &nproc);
  int icon = -2;
  Cblacs_get(-1, 0, &icon);
  Cblacs_gridinit(&icon, (char *) "r", gridrow, gridcol);

  test_svd(rank, nproc, icon);

  // exit
  MPI_Barrier(MPI_COMM_WORLD);
  Cblacs_exit(0);
}
