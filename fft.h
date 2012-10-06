#ifndef FFT_H
#define FFT_H

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <string>


using namespace std;

void inline bitreverse(vector<complex<double> > &sig) {
  unsigned int len = sig.size();
  unsigned int N = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(len))/log10(2.0)));
  unsigned int rev = 0;
// Processing Input Data
  for (unsigned int iter = 0; iter < N; ++iter)
    {
      if (rev > iter)
	{
// Replacing current values with reversed values

	  double tempr = real(sig[rev]);
	  double tempi = imag(sig[rev]);
	  complex<double> temp(tempr,tempi);
	  sig[rev] = sig[iter];
	  sig[iter] = temp;

	}
// Using filter "filt" such that the value of reverse changes with each iteration
      unsigned int filt = N;
      while (rev & (filt >>= 1)) {
	rev &= ~filt;
      }
      rev |= filt;
    }
  

}

void inline fft(vector<complex<double> > &data, int sign,unsigned int N){
  double pi = - 3.14159265358979;
  if ( sign == 1 || sign == -1) {
    pi = sign * pi;
  } else {
    cout << "Format fft(data, num), num = +1(fft) and num = -1 (Ifft)" << endl;
    exit(1);
  }
  unsigned int len = data.size();
  vector<complex<double> >::iterator it;
  it = data.end();
  if ( len != N) {
    unsigned int al = N - len;
    data.insert(it,al,complex<double>(0,0));
  }

  unsigned int K = (unsigned int) pow(2.0,ceil(log10(static_cast<double>(N))/log10(2.0)));
  vector<complex<double> >::iterator it1;
  it1 = data.end();
  if ( N < K) {
    unsigned int al = K - N;
    data.insert(it1,al,complex<double>(0,0));
    N = K;
  }

  bitreverse(data);

// radix2(data);
  for (unsigned int iter = 1; iter < N; iter <<=1)
    {
      const unsigned int step = iter << 1;

      const double theta = pi / double(iter);

      double wtemp = sin(theta * .5);
// Multipliers
      double wreal = -2 * wtemp * wtemp;
      double wimag = sin(theta);

// Factors
      double wr = 1.0;
      double wi = 0.0;
// Iteration through two loops

      for (unsigned int m = 0; m < iter; m++)
	{
// Iteration within m
	  for (unsigned int i = m; i < N; i += step)
	    {
// jth position
	      const unsigned int j = i + iter;

	      double tempr= wr * real(data[j]) - wi * imag(data[j]);
	      double tempi= wr * imag(data[j]) + wi * real(data[j]);

	      complex<double> temp(tempr,tempi);
	      data[j]= data[i]- temp;
	      data[i] += temp;

	    }
// Twiddle Factors updated
	  wtemp = wr;
	  wr += wr * wreal - wi * wimag;
	  wi += wi * wreal + wtemp * wimag ;
	}

    }

  if ( sign == -1) {
    double scale = 1.0/double(N);
    for (unsigned int i = 0; i < N; i++){
      data[i]*=scale;
    }
  }



// Place holder
  
}



void inline convfft(vector<double> &a, vector<double> &b, vector<double> &c) {
     unsigned int N = a.size() + b.size() - 1;
     vector<complex<double> > inp, filt;
     for (unsigned int i =0; i < a.size(); i++) {
         double temp = a[i];
         inp.push_back(complex<double>(temp,0));
     }
     for (unsigned int i =0; i < b.size(); i++) {
         double temp = b[i];
         filt.push_back(complex<double>(temp,0));
     }
     fft(inp,1,N);
     fft(filt,1,N);
     vector<complex<double> > temp;
     unsigned int K=inp.size();
     for (unsigned int i =0; i < K; i++){
         complex<double> mult = inp[i]*filt[i];
         temp.push_back(mult);
     }
     fft(temp, -1 , K);
     for (unsigned int i =0; i < N ; i++){
         double temp1 =real(temp[i]);
         c.push_back(temp1);
     }

}

#endif //FFT_H
