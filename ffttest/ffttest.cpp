#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <complex>
#include <vector>
#include <string>
#include "cycle.h"


using namespace std;

#define PI 3.14159265

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

void TSHkj(vector< complex<double> > &x, int q, int sgn, vector<double> &wrlong, vector<double> &wilong) {
  int n = x.size();
  int L = (int) pow(2.0, (double)q);
  int Ls = L / 2;
  int r = n / L;
  vector<double> wlr,wli;
  int t = (int) ceil(log10(static_cast<double>(n))/log10(2.0));
  int step = (int) pow((double)2, (double) t-q);
  for (int j = 0; j < Ls; j++) {
    wlr.push_back(wrlong[j*step]);
    wli.push_back(wilong[j*step]);
  }
  complex<double> tau;
  vector< complex<double> > y = x;
  

  for (int k = 0; k < r; k++) {
    for (int j = 0; j < Ls; j++) {
      double xlsr = real( y[(k+r)*Ls+j]);
      double xlsi = imag( y[(k+r)*Ls+j]);
      tau =complex<double>( wlr[j] * xlsr - wli[j] * xlsi,wlr[j] * xlsi + wli[j] * xlsr );
      x[k*L+j+Ls] = y[k*Ls+j] - tau;
      x[k*L+j] = y[k*Ls+j] + tau;
    }
      
  }
}

void AQkj(vector< complex<double> > &x, int q, int sgn, vector<double> &wrlong, vector<double> &wilong) {
  int n = x.size();
  int L = (int) pow(2.0, (double)q);
  int Ls = L / 2;
  int r = n / L;
  vector<double> wlr,wli;
  int t = (int) ceil(log10(static_cast<double>(n))/log10(2.0));
  int step = (int) pow((double)2, (double) t-q);
  for (int j = 0; j < Ls; j++) {
    wlr.push_back(wrlong[j*step]);
    wli.push_back(wilong[j*step]);
  }
  complex<double> tau;

  for (int k = 0; k < r; k++) {
    for (int j = 0; j < Ls; j++) {
      double xlsr = real( x[k*L+j+Ls]);
      double xlsi = imag( x[k*L+j+Ls]);
      tau =complex<double>( wlr[j] * xlsr - wli[j] * xlsi,wlr[j] * xlsi + wli[j] * xlsr );
      x[k*L+j+Ls] = x[k*L+j] - tau;
      x[k*L+j] = x[k*L+j] + tau;
    }
      
  }
}

void SHjk(vector< complex<double> > &x, int q, int sgn) {
  int n = x.size();
  int L = (int) pow(2.0, (double)q);
  int Ls = L / 2;
  int r = n / L;
  int rs = n / Ls;
  //vector<complex<double> > wl;
  double pi;
  if (sgn == 1 || sgn == -1) {
    pi = -1.0 * PI * sgn;
  } else {
    cout << "Please enter either 1(FFT) or -1(IFFT) for integer value sgn." << endl;
  }
  
  complex<double> tau;
  vector<complex<double> > y = x;

  double theta = 2*pi/L;
  double S = sin(theta);
  double C = cos(theta);
  double wlr = 1.0;
  double wli = 0.0;

  for (int j = 0; j < Ls; j++) {
    //wl.push_back(complex<double>(cos(2*PI*j/L),sn *sin(2*PI*j/L) ));
    //double theta = 2*PI*j/L;
    //double wlr = cos(theta);
    //double wli = sn * sin(theta);
    for (int k = 0; k < r; k++) {
      double xlsr = real( y[j*rs+k+r]);
      double xlsi = imag( y[j*rs+k+r]);
      tau =complex<double>( wlr * xlsr - wli * xlsi,wlr * xlsi + wli * xlsr ) ;
      x[(j+Ls)*r+k] = y[j*rs+k] - tau;
      x[j*r+k] = y[j*rs+k] + tau;
    }
    double temp = wlr;
    wlr = C * wlr - S * wli;
    wli = S * temp + C * wli;
      
  }
}

void GSjk(vector< complex<double> > &x, int q, int sgn) {
  int n = x.size();
  int L = (int) pow(2.0, (double)q);
  int Ls = L / 2;
  int r = n / L;
  //vector<complex<double> > wl;
  //double sn = (double) -1.0 * sgn;
  double pi;
  if (sgn == 1 || sgn == -1) {
    pi = -1.0 * PI * sgn;
  } else {
    cout << "Please enter either 1(FFT) or -1(IFFT) for integer value sgn." << endl;
  }
    

  complex<double> tau;
  double theta = 2*pi/L;
  double S = sin(theta);
  double C = cos(theta);
  double wlr = 1.0;
  double wli = 0.0;

  for (int j = 0; j < Ls; j++) {
    //wl.push_back(complex<double>(cos(2*PI*j/L),sn *sin(2*PI*j/L) ));
    //double theta = 2*PI*j/L;
    //double wlr = cos(theta);
    //double wli = sn * sin(theta);
    for (int k = 0; k < r; k++) {
      tau = x[k*L+j+Ls];
      double xlsr = real( x[k*L+j] - tau);
      double xlsi = imag( x[k*L+j] - tau);
      x[k*L+j+Ls] =complex<double>( wlr * xlsr - wli * xlsi,wlr * xlsi + wli * xlsr ) ;
      x[k*L+j] = x[k*L+j] + tau;
    }
    
    double temp = wlr;
    wlr = C * wlr - S * wli;
    wli = S * temp + C * wli;
    
      
  }
}

void AQjk(vector< complex<double> > &x, int q, int sgn) {
  int n = x.size();
  int L = (int) pow(2.0, (double)q);
  int Ls = L / 2;
  int r = n / L;
  //vector<complex<double> > wl;
  //double sn = (double) -1.0 * sgn;
  double pi;
  if (sgn == 1 || sgn == -1) {
    pi = -1.0 * PI * sgn;
  } else {
    cout << "Please enter either 1(FFT) or -1(IFFT) for integer value sgn." << endl;
  }
    

  complex<double> tau;
  double theta = 2*pi/L;
  double S = sin(theta);
  double C = cos(theta);
  double wlr = 1.0;
  double wli = 0.0;

  for (int j = 0; j < Ls; j++) {
    //wl.push_back(complex<double>(cos(2*PI*j/L),sn *sin(2*PI*j/L) ));
    //double theta = 2*PI*j/L;
    //double wlr = cos(theta);
    //double wli = sn * sin(theta);
    for (int k = 0; k < r; k++) {
      double xlsr = real( x[k*L+j+Ls]);
      double xlsi = imag( x[k*L+j+Ls]);
      tau =complex<double>( wlr * xlsr - wli * xlsi,wlr * xlsi + wli * xlsr ) ;
      x[k*L+j+Ls] = x[k*L+j] - tau;
      x[k*L+j] = x[k*L+j] + tau;
    }
    
    double temp = wlr;
    wlr = C * wlr - S * wli;
    wli = S * temp + C * wli;
    
      
  }
}

void fftct(vector<complex<double> > &data,int sgn, unsigned int N) {
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
  int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));

  for (int i=1; i < t+1; i++) {
    AQjk(data,i,sgn);
    
  }

}

void fft_sh(vector<complex<double> > &data,int sgn, unsigned int N) {
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

  int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));

  for (int i=1; i < t+1; i++) {
    SHjk(data,i,sgn);
    
  }

}

void fftct2(vector<complex<double> > &data,int sgn, unsigned int N) {
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
  int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));

  int L = (int) pow(2.0, (double)t);
  int Ls = L / 2;
  vector<double> wrlong,wilong;
  double sn = (double) -1.0 * sgn;
  for (int j = 0; j < Ls; j++) {
    wrlong.push_back(cos(2*PI*j/L));
    wilong.push_back(sn * sin(2*PI*j/L));
  }

  for (int i=1; i < t+1; i++) {
    AQkj(data,i,sgn,wrlong,wilong);
    
  }

}

void fft_tsh(vector<complex<double> > &data,int sgn, unsigned int N) {
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

  int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));

  int L = (int) pow(2.0, (double)t);
  int Ls = L / 2;
  vector<double> wrlong,wilong;
  double sn = (double) -1.0 * sgn;
  for (int j = 0; j < Ls; j++) {
    wrlong.push_back(cos(2*PI*j/L));
    wilong.push_back(sn * sin(2*PI*j/L));
  }

  for (int i=1; i < t+1; i++) {
    TSHkj(data,i,sgn,wrlong,wilong);
    
  }

}

void fftgs(vector<complex<double> > &data,int sgn, unsigned int N) {
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

  int t = (int) ceil(log10(static_cast<double>(N))/log10(2.0));

  for (int i=t; i > 0; i--) {
    GSjk(data,i,sgn);
    
  }

  bitreverse(data);

}


void convfft(vector<double> &a, vector<double> &b, vector<double> &c) {
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


int main() {
  int N = 64;
  vector<complex<double> > signal;
  for (int i =0; i < N; i++){
    signal.push_back(complex<double>((double)i, 0.0));
    cout << real(signal[i]) << " " << imag(signal[i]) << endl;
  }
  vector<complex<double> > sig1,sig2;
  sig1=signal;
  sig2=signal;
  
  double tdl, tct;
  ticks t0 = getticks();

  fftgs(sig1,1,N);
  ticks t1 = getticks();
  ticks t2 = getticks();
  fftct(sig2,1,N);
  ticks t3 = getticks();

  tdl = elapsed(t1,t0);
  tct = elapsed(t3,t2);
  
  for (int i =0; i < N; i++){
    cout << real(sig2[i])-real(sig1[i]) << " " << imag(sig2[i])-imag(sig1[i]) << endl;
  }
  //IFFT
  fftgs(sig1,-1,N);
  fftct(sig2,-1,N);
  cout << "IFFT - signal" << endl;
  for (int i =0; i < N; i++){
    cout << (real(sig2[i]) / N) << " " << (imag(sig2[i])/N) << endl;
    }

  cout << "tdl : " << tdl << "  " << "tct : " << tct << endl;
  
  return 0;
}
