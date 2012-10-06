#ifndef STATS_H
#define STATS_H

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <algorithm>
#include "fft.h"


using namespace std;

template <typename T>
T inline mean(vector<T> &vec) {
  T temp = (T) 0.0;
  int N = vec.size();
  for (int i=0; i < N ; i++) {
    temp+=vec[i];
    cout << temp << endl;
  }
 
  temp = (T) temp/N;
  return temp;
}

template <typename T>
int inline xcovar(vector<T> &x, vector<T> &y, vector<T> &z) {
  int lenx = x.size();
  int leny = y.size();
  int N;
  if (lenx >= leny){
    N = lenx;
    for (int j = 0; j < (lenx - leny); j++){
      y.push_back((T) 0.0);
    }
  } else {
    N = leny;
    for (int j = 0; j < (leny - lenx); j++) {
      x.push_back((T) 0.0);
    }
  }
  
  vector<double> a,b,c;
  T xm = mean(x);
  T ym = mean(y);

  for (int i = 0; i < lenx; i++)
    a.push_back((double) (x[i] - xm));
  for (int i = 0; i < leny; i++)
    b.push_back((double) (y[i] - ym));

  reverse(a.begin(),a.end());
  convfft(a,b,c);

  for (int i = N-1; i < (int) c.size();i++)
    z.push_back((T) c[i]/(double)N);

  return N;
}

template <typename T>
void xcovar(vector<T> &x, vector<T> &y, vector<T> &z, int lag) {
  vector<T> z2;
  int N;
  N = xcovar(x,y,z2);
  if (lag > N-1) {
    lag = N-1;
  }
  for (int i= 0; i < lag+1; i++) {
    z.push_back(z2[i]);
  }
}

template <typename T>
void inline xcorr(vector<T> &x, vector<T> &y, vector<T> &z) {
  xcovar(x,y,z);
  T temp = z[0];
  for (int i=0; i < (int) z.size(); i++) {
    z[i] = z[i]/temp;
  }
}

template <typename T>
void inline xcorr(vector<T> &x, vector<T> &y, vector<T> &z, int lag) {
  xcovar(x,y,z,lag);
  T temp = z[0];
  for (int i=0; i < (int) z.size(); i++) {
    z[i] = z[i]/temp;
  }
}

template <typename T>
void inline autocovar(vector<T> &x, vector<T> &z) {
  xcovar(x,x,z);
}

template <typename T>
void inline autocovar(vector<T> &x, vector<T> &z, int lag) {
  xcovar(x,x,z,lag);
}

template <typename T>
void inline autocorr(vector<T> &x, vector<T> &z) {
  xcorr(x,x,z);
}

template <typename T>
void inline autocorr(vector<T> &x, vector<T> &z, int lag) {
  xcorr(x,x,z,lag);
}
#endif //STATS_H
