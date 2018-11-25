#pragma once
#include <vector>
#include <iostream>
#include <complex>

using namespace std;

extern double GL_B,GL_W_A,GL_W_F;

class KetVector{
  vector<int> data;
  int energy_level;
  int size;
public:
  //создает вектор состоящий из нулей
  //н-р:|0000...0>
  KetVector(int size,int energy_level);
  //вычисляет сумму энергии
  //н-р:|2001> -> 3
  int getEnergysum();
  void Print();
  //вычисляет есть ли атомы
  //н-р:|2001> -> true
  //н-р:|2000> -> false
  bool haveNoAtoms();
  //находит следующий подходящий вектор
  //с такой же суммой энергии
  //н-р:|2001> -> |2010>
  //н-р:|2100> -> |3000>
  void next();
  //находит последний вектор
  //н-р: для случая с энергией == 3 и размером == 4
  //результат -> |3000>
  KetVector last();
  //прибавляет единицу к вектору,не меняя фотоны
  //н-р:|2010> -> |2011>
  //н-р:|2111> -> |2000>
  KetVector operator ++ (int);
  bool operator == (KetVector& x);
  bool operator != (KetVector& x);
  //возвращает i-ый элемент вектора,его можно поменять
  int& operator [] (int i);
  //возвращает результат в ячейке матрицы
  //н-р: для H4 и вектора размера == 4
  //н-р: i соотв-т |2011>
  //н-р: j соотв-т |1111>
  //результат в ячейке (i,j) == sqrt(2)*b;
  complex<double> operator * (KetVector y);
  //возращает в скольких атомах он различается с y
  int AtomDifference(KetVector y);
//end of class
};
