#include "DistributedMatrix.h"

DistributedMatrix::DistributedMatrix(int rows,int cols,int lrank,int lnproc,int lblock,int licon) {
  rank = lrank;
  nproc = lnproc;
  icon = licon;
  blocksize = lblock;
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

int* DistributedMatrix::descinit() {
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

void DistributedMatrix::Create() {
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

void DistributedMatrix::Fill(string name) {
  MPI_File thefile;
  MPI_Offset filesize;
  MPI_File_open(MPI_COMM_WORLD,name.c_str(),MPI_MODE_RDONLY,MPI_INFO_NULL,&thefile);
  MPI_File_get_size(thefile,&filesize);
  filesize = filesize/sizeof(double);
  const int matrix = rows*cols*2,qvector = rows*2,number = 1;
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
      double* tmp_buf = new double[qvector];
      MPI_File_read(thefile,tmp_buf,qvector,MPI_DOUBLE,MPI_STATUS_IGNORE);
      for (int i = 0; i < myProcRows; i++) {
        complex<double> ith(tmp_buf[(i+myProcRow*blocksize)*2],tmp_buf[(i+myProcRow*blocksize)*2+1]);
        for (int j = 0; j < myProcCols; j++) {
          complex<double> x(tmp_buf[i*2],tmp_buf[i*2+1]);
          complex<double> y(tmp_buf[j*2],tmp_buf[j*2+1]);
          // data[i*myProcCols + j] = x*y;

          complex<double> jth(tmp_buf[(j+myProcCol*blocksize)*2],tmp_buf[(j+myProcCol*blocksize)*2+1]);
          data[i*myProcCols + j] = ith*jth;
        }
      }
      delete(tmp_buf);
    }
    else{
      cout<<filesize<<endl;
      cout<<"ERROR::we are not prepared"<<endl;
      /*
      double* tmp_buf = new double[number];
      MPI_File_read(thefile,tmp_buf,number,MPI_DOUBLE,MPI_STATUS_IGNORE);
      for (int i = 0; i < myProcRows; i++) {
        for (int j = 0; j < myProcCols; j++) {
          if(myProcCol == ind)
            data[i*myProcCols + j] = 0;
        }
      }
      */
    }
  }
}

void DistributedMatrix::Fill(complex<double> c) {
  for (int i = 0; i < myProcRows; i++) {
    for (int j = 0; j < myProcCols; j++) {
      data[i*myProcCols + j] = c;
    }
  }
}

void DistributedMatrix::Fill(vector<KetVector> v) {
  for (int i = 0; i < myProcRows; i++) {
    for (int j = 0; j < myProcCols; j++) {
      int v_i = myProcRowsOffset + i;
      int v_j = myProcColsOffset + j;
      complex<double> c = v[v_i]*v[v_j];
      data[i*myProcCols + j] = c;
    }
  }
}

/*
 * заполняет диагональ матрицы значениями из double*
 */
void DistributedMatrix::Fill(double* x){
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

void DistributedMatrix::Print() {
  for(int i = 0; i < procrows; ++i) {
    for(int j = 0; j < proccols; ++j) {
      MPI_Barrier(MPI_COMM_WORLD);
      if(i*proccols + j == rank) {
        cout<<"coordinates of processes grid:"<<endl;
        cout<<"("<<myProcRow<<","<<myProcCol<<")"<<endl;
        for(int x = 0; x < myProcRows; ++x) {
          for(int y = 0; y < myProcCols; ++y) {
            cout<< "["<<data[x*myProcCols + y].real()<<"]";
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

void DistributedMatrix::Print(string key) {
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

complex<double>* DistributedMatrix::getData() {
  return data;
}

int* DistributedMatrix::getDesc() {
  return desc;
}

/*
 * result: C=C+this*B
 * C изменится,
 * C меняется во время процедуры, так что не использовать ее в качестве this или B
 */
void DistributedMatrix::Multiply(DistributedMatrix& B, DistributedMatrix& C) {
  complex<double> cdOne(1.0,0.0);
  complex<double> cdZero(0.0, 0.0);
  int ione = 1;

  pzgemm_((char*)"N", (char*)"N", &rows, &rows, &rows,
          &cdOne, data, &ione, &ione, desc,
          B.getData(), &ione, &ione, B.getDesc(),
          &cdZero, C.getData(), &ione, &ione, C.getDesc());
}

void DistributedMatrix::Multiply(DistributedMatrix& A,DistributedMatrix& B,DistributedMatrix& C){
  DistributedMatrix tmp1(rows,cols,rank,nproc,blocksize,icon);
  DistributedMatrix tmp2(rows,cols,rank,nproc,blocksize,icon);
  B.Multiply(C,tmp1);
  A.Multiply(tmp1,tmp2);
  *this = tmp2;
}

void DistributedMatrix::conjTrans(DistributedMatrix& B) {
  complex<double> cdOne(1.0,0.0);
  complex<double> cdZero(0.0, 0.0);
  int ione = 1;
  DistributedMatrix I(rows,cols,rank,nproc,blocksize,icon);
  double* ones = new double[rows];
  for(int i = 0;i < rows;++i) {
    ones[i] = 1;
  }
  I.Fill(ones);
  pzgemm_((char*)"C", (char*)"N", &rows, &rows, &rows,
          &cdOne, data, &ione, &ione, desc,
          I.getData(), &ione, &ione, I.getDesc(),
          &cdZero, B.getData(), &ione, &ione, B.getDesc());
}

void DistributedMatrix::CreateSVD(DistributedMatrix& middle, DistributedMatrix& u, DistributedMatrix& vt) {
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

void DistributedMatrix::Diagonalize(DistributedMatrix& middle, DistributedMatrix& u, DistributedMatrix& vt) {
  double* s = new double [rows];

  int iOne = 1;
  complex<double> *work = new complex<double>[1];
  int info = 0, lwork = -1, lrwork = -1, liwork = -1;
  double *rwork = new double[1];
  int *iwork = new int[1];

  pzheevd_((char*) "V", (char*) "U",
    &rows, data, &iOne, &iOne, desc, s,
    u.getData(), &iOne, &iOne, u.getDesc(),
    work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

  lwork = (int) work[0].real();
  delete[] work;
  work = new complex<double>[lwork];

  lrwork = (int) rwork[0];
  delete[] rwork;
  rwork = new double[lrwork];

  liwork = (int) iwork[0];
  delete[] iwork;
  iwork = new int[liwork];

  pzheevd_((char*) "V", (char*) "U",
    &rows, data, &iOne, &iOne, desc, s,
    u.getData(), &iOne, &iOne, u.getDesc(),
    work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

  u.conjTrans(vt);
  middle.Fill(s);
  delete(s);
}

void DistributedMatrix::Exp(double dt){
  if(myProcRow != myProcCol)
    return;
  for(int i = 0;i < myProcRows;++i){
    data[i*myProcCols + i] = exp(complex<double>(0,-1) * data[i*myProcCols + i] * dt);
  }
}

DistributedMatrix DistributedMatrix::operator = (DistributedMatrix& t){
  for(int i = 0;i < myProcRows;++i) {
    for(int j = 0;j < myProcCols;++j) {
      data[i*myProcCols + j] = t.getData()[i*myProcCols + j];
    }
  }
  return *this;
}
