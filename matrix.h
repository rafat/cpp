#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

template <typename T>
class mat {
  int mrows;
  int mcols;
  vector<T> mdata;

public:
  mat() {

  }
  mat(int r, int c) {
    mrows = r;
    mcols = c;
    mdata.resize(r*c);
    
}
  mat(int r, int c, vector<T> &datval) {
    if ( r*c == (int) datval.size()) {
      mrows = r;
      mcols = c;
      mdata = datval;
    } else {
      cout << "Matrix Dimensions and length of Data Vector Don't Match" << endl;
    }
  }

  T& operator()(int i,int j){
    return mdata[i * mcols + j];
}

  T operator()(int i,int j) const{
    return mdata[i * mcols + j];
}

  void zeros() {
    mdata.clear();
    mdata.insert(mdata.begin(),mrows*mcols, (T) 0.0);

  }

  void ones() {
    mdata.clear();
    mdata.insert(mdata.begin(),mrows*mcols, (T) 1.0);

  }

  int rows() {
    return mrows;   
  }

  int cols() {
    return mcols;
  }

  void disp() {
    for (int i=0; i < mrows; i++) {
      cout << "Row" << i << " : ";
      for (int j=0; j < mcols; j++) {
	cout << mdata[i*mcols+j] << " ";
      }
      cout << endl;
    }
  }

  void setMatrix(int r, int c, vector<T> &datval) {
    mrows = r;
    mcols = c;
    mdata = datval;
  }

  void getData(vector<T> &send_data) {
    send_data = mdata;
  }

  mat operator+ (mat &B) {
    mat<T> temp(mrows,mcols);
    if ( ( mrows == B.rows() ) && (mcols == B.cols())) {
      for (int i = 0; i < (int) B.rows();i++) {
	for (int j = 0; j < (int) B.cols();j++) {
	  temp(i,j) = mdata[i * mcols + j] + B(i,j);
	}
      }
    }
    else {
      cout << " Matrix Dimensions Must Agree" << endl;
    }
    return temp;
  }

  mat operator- (mat &B) {
    mat<T> temp(mrows,mcols);
    if ( ( mrows == B.rows() ) && (mcols == B.cols())) {
      for (int i = 0; i < (int) B.rows();i++) {
	for (int j = 0; j < (int) B.cols();j++) {
	  temp(i,j) = mdata[i * mcols + j] - B(i,j);
	}
      }
    }
    else {
      cout << " Matrix Dimensions Must Agree" << endl;
    }
    return temp;
  }

  mat operator* (mat &B) {
    mat<T> temp(mrows,B.cols());
    if ( mcols == B.rows()) {
      for (int i = 0; i < (int) mrows;i++) {
	for (int j = 0; j < B.cols(); j++) {
	  T temp2 =(T) 0.0;
	  for (int k = 0; k < B.rows();k++) {
	    temp2 = temp2 + mdata[i*mcols+k]* B(k,j);
	  }
	  temp(i,j) = temp2;
	}
      }
    }
    else {
      cout << " Matrix Dimensions Must Agree" << endl;
    }
    return temp;
  }

  void eye(int N) {
    mrows = N;
    mcols = N;
    zeros();
    for (int i = 0; i < N; i++) {
      mdata[i*mcols+i] = (T) 1.0;
    }
  }


  virtual ~mat() {
  }

};


#endif //DATA_H
