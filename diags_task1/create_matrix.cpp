#include <cstdio>
#include <string>
#include <cstdlib>
#include <complex>

using namespace std;

void function(int i,int j,double& real,double& img){
  complex<double> res;
  res = 6*i + j;
  real = res.real();
  img = res.imag();
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
