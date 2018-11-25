#include <cstdio>
#include <complex>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "mpi.h"
#include "scalapack.h"
#include "DistributedMatrix.h"
#include "KetVector.h"

using namespace std;

double GL_W_F = 1;
double GL_W_A = 2;
double GL_B = 3;


void Init_grid(int& icon,int rank,int nproc,int& blocksize,int dim){
  int gridrow,gridcol;
  int divA;
  for (divA = (int) sqrt(nproc); (nproc % divA != 0); divA--);
  int divB = nproc / divA;
  gridrow = divB;
  gridcol = divA;
  Cblacs_gridinit(&icon, (char *) "r", gridrow, gridcol);
  blocksize = ceil(1.0*dim/gridcol);
  if(rank == 0){
    cout<<"----------------"<<endl;
    cout<<"grid of processes: "<<gridrow<<"x"<<gridcol<<endl;
    cout<<"blocksize: "<<blocksize<<endl;
    cout<<"----------------"<<endl;
  }
}

//создает последовательность векторов(те,что мы пишем сверху и слева матрицы)
vector<KetVector> genVectors(int size,int energy_level){
  KetVector x(size,energy_level);
  KetVector last = x.last();
  vector<KetVector> result;
  for(x.next();x != last;x.next()){
    result.push_back(x);
  }
  result.push_back(x);
  return result;
}

//выводит последовательно H1,H2,...,Hn
//тебе надо как то написать вывод большой матрицы.
//я изменил вывод в DistributedMatrix::Print(),
//чтобы выводил только real часть
void test_task2(int rank,int nproc){
  int size = 4;
  int* icon = new int[size];
  int blocksize;
  for(int energy_level = 1;energy_level <= size;++energy_level){
    vector<KetVector> v = genVectors(size,energy_level);
    int dim = v.size();
    Cblacs_get(-1, 0, &(icon[energy_level - 1]));
    Init_grid(icon[energy_level - 1],rank,nproc,blocksize,dim);
    DistributedMatrix H(dim,dim,rank,nproc,blocksize,icon[energy_level - 1]);
    H.Fill(v);
    H.Print();
  }
}

//надо сделать ввод данных ч/з аргументы кс
//w_a,w_c,b(те,что внутри матрицы) определены как глобальные переменные(вверху main.cpp)
int main(int argc,char** argv){
  // init
  int rank,nproc,blocksize,dim;

  MPI_Init(&argc, &argv);
  Cblacs_pinfo(&rank, &nproc);

  test_task2(rank,nproc);

  // exit
  MPI_Barrier(MPI_COMM_WORLD);
  Cblacs_exit(0);
}
