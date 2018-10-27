#include <cstdio>
#include <string>
#include <cstdlib>
#include <complex>

using namespace std;

/*
 * 1 1 1 2 2 2
 * 1 1 1 2 2 2
 * 1 1 1 2 2 2
 * 3 3 3 4 4 4
 * 3 3 3 4 4 4
 * 3 3 3 4 4 4
 */
void block(int i, int j, double& real, double& img) {
  complex<double> res;
  if((i>=3) && (j<3))
    real = 3;
  else{
    if((i>=3) && (j>=3))
      real = 4;
    if((i<3) && (j>=3))
      real = 2;
    if((i<3) && (j<3))
      real = 1;
  }
  /*
  if(j>=3)
    real = 2;
  //res = 6*i + j;
  */
  img = 0;
}

/*
 * 0  1  2  3  4  5
 * 6  7  8  9  10 11
 * 12 13 14 15 16 17
 * 18 19 20 21 22 23
 * 24 25 26 27 28 29
 * 30 31 32 33 34 35
 */
void line(int i, int j, double& real, double& img) {
  complex<double> res;
  res = 6*i + j;
  real = res.real();
  img = res.imag();
}


void function(int i,int j,double& real,double& img){
  line(i, j, real, img);
}

void create_file(size_t N,size_t M,char* str){
  FILE* f;
  f = fopen(str,"wb");
	for(size_t i = 0;i < N;++i)
  	for(size_t j = 0;j < M;++j){
      double real,img;
    	function(i,j,real,img);
    	fwrite(&real,sizeof(real),1,f);
    	fwrite(&real,sizeof(img),1,f);
	   }
  fclose(f);
}

int main(int argc,char** argv){
	size_t x,y;
	x = atoi(argv[1]);
	y = atoi(argv[2]);
	create_file(x,y,argv[3]);
	return 0;
}
