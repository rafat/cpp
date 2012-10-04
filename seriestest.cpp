#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "series.h"

using namespace std;

int main() {
  int row = 10;
  int col = 10;
  mat<double> matrix(row,col);
  for (int i=0; i < row; i++) {
    for (int j=0; j < col; j++) {
      matrix(i,j) =(double) i*col +j;
    }
  }
  series<double> dseries(matrix);
  
  matrix.disp();
  vector<double> datavec;
  dseries.getData(datavec);
  for (int i = 0; i < (int) datavec.size(); i++)
    cout << datavec[i] << " ";
  return 0;

}
