#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "matrix.h"

using namespace std;

int main() {
  int row = 4;
  int col = 5;
  int row2 = 5;
  int col2 = 6;
  mat<double> matrix(row,col);
  mat<double> matrix2(row,col);
  mat<double> matrix3(row2,col2);
  mat<double> matsum;
  for (int i=0; i < row; i++) {
    for (int j=0; j < col; j++) {
      matrix(i,j) = i*col +j;
    }
  }
  cout << matrix(4,4) << endl;
  matrix2 = matrix;
  // matrix.ones();
  //cout << matrix(4,4) << endl;
  matrix.disp();
  vector<double> getvec;
  matrix.getData(getvec);
  for (int i = 0; i < (int) getvec.size(); i++)
    cout << getvec[i] << " " << endl;
  matsum = matrix + matrix2;
  matsum.disp();

  for (int i=0; i < row2; i++) {
    for (int j=0; j < col2; j++) {
      matrix3(i,j) = i*col2 +j;
    }
  }
  matrix3.disp();

  mat<double> matmult;
  matmult = matrix * matrix3;
  matmult.disp();
  mat<double> etest;
  etest.eye(4);
  etest.disp();

  vector<double> vec1;
  vec1.push_back(1.0);
  vec1.push_back(7.0);
  vec1.push_back(6.2);

  etest.diag(vec1);
  etest.disp();
  
  return 0;

}
