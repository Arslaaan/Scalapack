#pragma once

#include "mpi.h"
#include "scalapack.h"
#include "KetVector.h"

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
  DistributedMatrix(int rows,int cols,int lrank,int lnproc,int lblock,int licon);
  int* descinit();
  void Create() ;
  void Fill(string name);
  void Fill(complex<double> c);
  void Fill(vector<KetVector> v);
  /*
   * заполняет диагональ матрицы значениями из double*
   */
  void Fill(double* x);
  void Print();
  void Print(string key);
  complex<double>* getData();
  int* getDesc();
  /*
   * result: C=C+this*B
   * C изменится,
   * C меняется во время процедуры, так что не использовать ее в качестве this или B
   */
  void Multiply(DistributedMatrix& B, DistributedMatrix& C);
  void Multiply(DistributedMatrix& A,DistributedMatrix& B,DistributedMatrix& C);
  void conjTrans(DistributedMatrix& B);
  void CreateSVD(DistributedMatrix& middle, DistributedMatrix& u, DistributedMatrix& vt);
  void Diagonalize(DistributedMatrix& middle, DistributedMatrix& u, DistributedMatrix& vt);
  void Exp(double dt);
  DistributedMatrix operator = (DistributedMatrix& t);
};
