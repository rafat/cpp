#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "matrix.h"

using namespace std;

int main() {
  int row = 10;
  int col = 10;
  mat<double> matrix(row,col);
  for (int i=0; i < row; i++) {
    for (int j=0; j < col; j++) {
      matrix(i,j) = i*col +j;
    }
  }
  cout << matrix(4,4) << endl;
  // matrix.ones();
  //cout << matrix(4,4) << endl;
  matrix.disp();
  return 0;

}
