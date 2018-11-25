#include "KetVector.h"

KetVector::KetVector(int size,int energy_level){
  this->energy_level = energy_level;
  this->size = size;
  for (size_t i = 0; i < size; i++) {
    data.push_back(0);
  }
}

int KetVector::getEnergysum(){
  int sum = 0;
  for (size_t i = 0; i < size; i++) {
    sum += data[i];
  }
  return sum;
}

void KetVector::Print(){
  cout<<"|";
  for (size_t i = 0; i < size; i++) {
    cout<<data[i];
  }
  cout<<">"<<endl;
}

bool KetVector::haveNoAtoms(){
  for (size_t i = 1; i < size; i++) {
    if(data[i] == 1)
      return false;
  }
  return true;
}

void KetVector::next(){
  KetVector x = *this;
  int es;
  x++;
  if(x.haveNoAtoms())
    x[0] += 1;
  while((es = x.getEnergysum()) != energy_level){
    x++;
    if(x.haveNoAtoms())
      x[0] += 1;
  }
  *this = x;
}

KetVector KetVector::last(){
  KetVector last(size,energy_level);
  last[0] = energy_level;
  return last;
}

KetVector KetVector::operator ++ (int){
  int i = size - 1;
  while(i > 0){
    if(data[i] == 1)
      data[i] = 0;
    else{
      data[i] = 1;
      break;
    }
    --i;
  }
  return *this;
}

bool KetVector::operator == (KetVector& x){
  for (size_t i = 0; i < size; i++) {
    if(data[i] != x[i])
      return false;
  }
  return true;
}

bool KetVector::operator != (KetVector& x){
  return !(*this == x);
}

int& KetVector::operator [] (int i){
  return data[i];
}

complex<double> KetVector::operator * (KetVector y){
  if(*this == y){
    double x = sqrt(y[0])*GL_W_F + sqrt(y.getEnergysum() - y[0])*GL_W_A;
    return complex<double> (x,0.0);
  }
  else{
    if((abs((*this)[0] - y[0]) == 1) && (this->AtomDifference(y) == 1)){
      double x = sqrt(max((*this)[0],y[0]))*GL_B;
      return complex<double> (x,0.0);
    }
    else
      return complex<double> (0.0,0.0);
  }
  return complex<double> (-1.0,0.0);
}

int KetVector::AtomDifference(KetVector y){
  int sum = 0;
  for (size_t i = 1; i < size; i++) {
    sum += abs(data[i] - y[i]);
  }
  return sum;
}
