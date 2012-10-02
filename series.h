#ifndef SERIES_H
#define SERIES_H

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "matrix.h"

using namespace std;

template <Typename T>
class series {
  vector<T> data;
  vector<string> index;
  vector<int> pseudo;
  int dimension = 1;

public:
  series(vector<T> &signal) {

    int lensig = signal.size();
    for (int i = 0; i < lensig; i++) {
      pseudo.push_back(i);

    }
    
    data = signal;
  }

  series(vector<T> &signal, vector<string> &ind) {
    int lensig = signal.size();
    for (int i = 0; i < lensig; i++) {
      pseudo.push_back(i);

    }
    
    data = signal;
    index = ind;


  }

  series(mat<T> &signal) {
    dimension = signal.rows();
    vector<T> sigvec;
    data = signal.getData(sigvec);
    for (int i = 0; i < signal.cols(); i++) {
      pseudo.push_back(i);
    }
   
  }

  series(mat<T> &signal, vector<string> &ind) {
    dimension = signal.rows();
    vector<T> sigvec;
    data = signal.getData(sigvec);
    for (int i = 0; i < signal.cols(); i++) {
      pseudo.push_back(i);
    }
    index = ind;
   
  }

  virtual ~series() {}

};


#endif //SERIES_H
