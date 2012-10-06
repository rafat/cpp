#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "stats.h"

using namespace std;

int main() {
  int row = 11;
  vector<double> inp;
  for (int i=0; i < row; i++)
    inp.push_back((double) i);

  cout << mean(inp) << endl;
  vector<double> oup;
  autocorr(inp,oup,5);
  for (int i=0; i < (int) oup.size(); i++) {
    cout << oup[i] << " ";
  }
  cout << endl;
  return 0;

}
