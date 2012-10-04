#ifndef SERIES_H
#define SERIES_H

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "matrix.h"

using namespace std;

template <typename T>
class series {
  vector<T> data;
  vector<string> index;
  vector<int> pseudo;
  int dimension;
  mat<T> matrix;

public:
  series(vector<T> &signal) {

    int lensig = signal.size();
    for (int i = 0; i < lensig; i++) {
      pseudo.push_back(i);

    }
    
    data = signal;
    matrix.setMatrix(1,lensig,data);
    dimension = 1;
  }

  series(vector<T> &signal, vector<string> &ind) {
    int lensig = signal.size();
    for (int i = 0; i < lensig; i++) {
      pseudo.push_back(i);

    }
    
    data = signal;
    index = ind;
    matrix.setMatrix(1,lensig,data);
    dimension = 1;
  }

  series(mat<T> &signal) {
    dimension = signal.rows();
    vector<T> sigvec;
    signal.getData(sigvec);
    data = sigvec;
    for (int i = 0; i < signal.cols(); i++) {
      pseudo.push_back(i);
    }
    matrix = signal;
   
  }

  series(mat<T> &signal, vector<string> &ind) {
    dimension = signal.rows();
    vector<T> sigvec;
    signal.getData(sigvec);
    data = sigvec;
    for (int i = 0; i < signal.cols(); i++) {
      pseudo.push_back(i);
    }
    index = ind;
    matrix = signal;
   
  }

  void getData(vector<T> &datavec) {
    datavec = data;
  }

  virtual ~series() {}

};


#endif //SERIES_H
