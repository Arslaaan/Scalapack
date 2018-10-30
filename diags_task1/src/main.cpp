#include <cstdio>
#include "mpi.h"
#include <complex>
#include <iostream>
#include "scalapack.h"
#include <iomanip>
#include <string>
#include "DistributedMatrix.cpp"

using namespace std;

void test_svd(int rank, int nproc, int icon,double dt,int n,int dim,int blocksize) {
  DistributedMatrix H(dim,dim,rank,nproc,blocksize,icon);
  DistributedMatrix rho(dim,dim,rank,nproc,blocksize,icon);
  DistributedMatrix tmp(dim,dim,rank,nproc,blocksize,icon);
  DistributedMatrix left(dim,dim,rank,nproc,blocksize,icon);
  DistributedMatrix right(dim,dim,rank,nproc,blocksize,icon);
  DistributedMatrix middle(dim,dim,rank,nproc,blocksize,icon);
  DistributedMatrix unitary(dim,dim,rank,nproc,blocksize,icon);
  DistributedMatrix transposed(dim,dim,rank,nproc,blocksize,icon);
  H.Fill("hamilton");
  rho.Fill("rho");
  left.Fill(static_cast<double>(0));
  middle.Fill(static_cast<double>(0));
  right.Fill(static_cast<double>(0));
  tmp.Fill(static_cast<double>(0));

  H.Diagonalize(middle,left,right);// H == left*middle*right

  middle.Exp(dt);
  unitary.Multiply(left,middle,right);//unitary = left*middle*right;

  unitary.conjTrans(transposed);


  for(int i = 0;i < n;++i){
    rho.Print("diagonal");
    tmp.Multiply(unitary,rho,transposed);
    rho = tmp;
  }
}


int main(int argc,char** argv){
  // init
  int rank,nproc;
  int gridrow = 2,gridcol = 2;
  double dt = atof(argv[1]);
  int n = atoi(argv[2]);
  int dim = atoi(argv[3]);
  MPI_Init(&argc, &argv);
  Cblacs_pinfo(&rank, &nproc);
  int icon = -2;
  Cblacs_get(-1, 0, &icon);
  int divA;
  for (divA = (int) sqrt(nproc); (nproc % divA != 0); divA--);
  int divB = nproc / divA;
  gridrow = divB;
  gridcol = divA;
  Cblacs_gridinit(&icon, (char *) "r", gridrow, gridcol);
  int blocksize = ceil(1.0*dim/gridcol);
  if(rank == 0){
    cout<<"----------------"<<endl;
    cout<<"grid of processes: "<<gridrow<<"x"<<gridcol<<endl;
    cout<<"blocksize: "<<blocksize<<endl;
    cout<<"----------------"<<endl;
  }
  test_svd(rank, nproc, icon,dt,n,dim,blocksize);

  // exit
  MPI_Barrier(MPI_COMM_WORLD);
  Cblacs_exit(0);
}
