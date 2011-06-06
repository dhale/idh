/****************************************************************************
Copyright (c) 2011, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
*****************************************************************************
Benchmark nested parallel loops with OpenMP.
Compile with: g++ -O3 -fopenmp ParallelNest.cpp.
@author Dave Hale, Colorado School of Mines
@version 2011.06.04
****************************************************************************/

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

/////////////////////////////////////////////////////////////////////////////
// stopwatch

#ifdef _WIN32
#include <ctime>
#else
#include <sys/time.h>
#endif

double const maxtime = 3.0;
int const maxtest = 3;

class Stopwatch {
public:
  Stopwatch() : running_(false), time_(0.0) {
  }
  void start() {
    if (!running_) {
      running_ = true;
      start_ = timeInSeconds();
    }
  }
  void stop() {
    if (running_) {
      time_ = timeInSeconds()-start_;
      running_ = false;
    }
  }
  void reset() {
    stop();
    time_ = 0.0;
  }
  void restart() {
    reset();
    start();
  }
  double time() const {
    if (running_) {
      return time_+timeInSeconds()-start_;
    } else {
      return time_;
    }
  }
private:
  bool running_;
  double start_;
  double time_;
#ifdef _WIN32
  double timeInSeconds() const {
    return (double)std::clock()/(double)CLOCKS_PER_SEC;
  }
#else
  double timeInSeconds() const {
    struct timeval t;
    gettimeofday(&t,0);
    return t.tv_sec+0.000001*t.tv_usec;
  }
#endif
};

/////////////////////////////////////////////////////////////////////////////
// utilities

typedef vector<float> vfloat1;
typedef vector<vfloat1> vfloat2;
typedef vector<vfloat2> vfloat3;

vfloat3 fillvfloat3(float f, int n1, int n2, int n3) {
  vfloat3 a(n3,vfloat2(n2,vfloat1(n1)));
  for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        a[i3][i2][i1] = f;
      }
    }
  }
  return a;
}

void assertEquals(const vfloat2& a, const vfloat2& b) {
  int n1 = a[0].size();
  int n2 = a.size();
  for (int i2=0; i2<n2; ++i2)
    for (int i1=0; i1<n1; ++i1)
      assert(a[i2][i1]==b[i2][i1]);
}
void assertEquals(const vfloat3& a, const vfloat3& b) {
  int n3 = a.size();
  for (int i3=0; i3<n3; ++i3)
    assertEquals(a[i3],b[i3]);
}
void assertEquals(int m, int n, float** a, float** b) {
  for (int i=0; i<m; ++i) {
    for (int j=0; j<n; ++j) {
      assert(a[i][j]==b[i][j]);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
// second-order linear recursion

void solr(
  float a1, float a2, float b0, float b1, float b2, 
  const vfloat1& x, vfloat1& y) {
  int n = x.size();
  float yim2 = 0.0f;
  float yim1 = 0.0f;
  float xim2 = 0.0f;
  float xim1 = 0.0f;
  for (int i=0; i<n; ++i) {
    float xi = x[i];
    float yi = b0*xi+b1*xim1+b2*xim2-a1*yim1-a2*yim2;
    y[i] = yi;
    yim2 = yim1;
    yim1 = yi;
    xim2 = xim1;
    xim1 = xi;
  }
}
void solrS(
  float a1, float a2, float b0, float b1, float b2,
  const vfloat2& x, vfloat2& y) {
  int n = x.size();
  for (int i=0; i<n; ++i)
    solr(a1,a2,b0,b1,b2,x[i],y[i]);
}
void solrP(
  float a1, float a2, float b0, float b1, float b2,
  const vfloat2& x, vfloat2& y) {
  int n = x.size();
  #pragma omp parallel for schedule(dynamic)
  for (int i=0; i<n; ++i)
    solr(a1,a2,b0,b1,b2,x[i],y[i]);
}
void solrS(
  float a1, float a2, float b0, float b1, float b2,
  const vfloat3& x, vfloat3& y) {
  int n = x.size();
  for (int i=0; i<n; ++i)
    solrS(a1,a2,b0,b1,b2,x[i],y[i]);
}
void solrP(
  float a1, float a2, float b0, float b1, float b2,
  const vfloat3& x, vfloat3& y) {
  int n = x.size();
  #pragma omp parallel for schedule(dynamic)
  for (int i=0; i<n; ++i)
    solrP(a1,a2,b0,b1,b2,x[i],y[i]);
}
void benchSolr(int n1, int n2, int n3) {
  cout<<"Array solr: n1="<<n1<<" n2="<<n2<<" n3="<<n3<<endl;
  int niter;
  double mflop2 = 9.0e-6*n1*n2;
  double mflop3 = 9.0e-6*n1*n2*n3;
  Stopwatch sw;
  float a1 = -1.8f;
  float a2 = 0.81f;
  float b0 = 2.0f;
  float b1 = -3.2f;
  float b2 = 1.28f;
  vfloat3 x = fillvfloat3(1.0f,n1,n2,n3);
  vfloat3 ys = fillvfloat3(0.0f,n1,n2,n3);
  vfloat3 yp = fillvfloat3(0.0f,n1,n2,n3);
  for (int ntest=0; ntest<maxtest; ++ntest) {
    sw.restart();
    for (niter=0; sw.time()<maxtime; ++niter)
      solrS(a1,a2,b0,b1,b2,x[niter%n3],ys[niter%n3]);
    sw.stop();
    cout<<"2D S: rate = "<<(int)(niter*mflop2/sw.time())<<endl;
    sw.restart();
    for (niter=0; sw.time()<maxtime; ++niter)
      solrP(a1,a2,b0,b1,b2,x[niter%n3],yp[niter%n3]);
    sw.stop();
    cout<<"2D P: rate = "<<(int)(niter*mflop2/sw.time())<<endl;
    assertEquals(ys[0],yp[0]);
    sw.restart();
    for (niter=0; sw.time()<maxtime; ++niter)
      solrS(a1,a2,b0,b1,b2,x,ys);
    sw.stop();
    cout<<"3D S: rate = "<<(int)(niter*mflop3/sw.time())<<endl;
    sw.restart();
    for (niter=0; sw.time()<maxtime; ++niter)
      solrP(a1,a2,b0,b1,b2,x,yp);
    sw.stop();
    cout<<"3D P: rate = "<<(int)(niter*mflop3/sw.time())<<endl;
    assertEquals(ys,yp);
  }
}

int main(int argc, char** argv) {
  int ns[5][3] = {
    {1000, 50000,     5},
    {1000,  5000,    50},
    {1000,   500,   500},
    {1000,    50,  5000},
    {1000,     5, 50000}};
  for (int i=0; i<5; i+=4) {
    int n1 = ns[i][0], n2 = ns[i][1], n3 = ns[i][2];
    benchSolr(n1,n2,n3);
  }
}
