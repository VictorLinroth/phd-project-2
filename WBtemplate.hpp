#ifndef WB_TEMPLATE
#define WB_TEMPLATE

#include <cmath>
#include <vector>
#include <memory>

template<typename T, typename F>
class WBcoeff {
private:
  std::vector<std::vector<T>> w;
  std::vector<std::vector<T>> A;
public:
  WBcoeff() : w(), A() {}
  WBcoeff(std::size_t N) : w(), A() {
    std::size_t m, n;
    find(m,n,N);
    extend_w(m,N);
  }
  void find(std::size_t &m,std::size_t &n,std::size_t N) {
    n = 0;
    while ((N%2)!=0) {
      N = (N-1)/2;
      ++n;
    }
    m = N/2;
  }

  void extend_w(std::size_t m, std::size_t N);
  void alloc_A(std::size_t m, std::size_t n) {
    if (A.size() < m+1) 
      A.resize(m+1);
    if (A[m].size() < n+1)
      A[m].resize(n+1,0.0);
  }
  void extend(std::size_t N) {
    std::size_t m, n;
    find(m,n,N);
    extend_w(m,N);
  }
  T getw(std::size_t m, std::size_t n) const {
    return w[m][n];
  }
  typename std::vector<T>::const_iterator cbeginw(std::size_t m) const {
    return w[m].cbegin();
  }
  std::size_t sizew(std::size_t m) {
    return w[m].size();
  }
  T getA(std::size_t m, std::size_t n) const {
    return A[m][n];
  }
  void setA(std::size_t m, std::size_t n, T val) {
    A[m][n] = val;
  }
};

template<typename T,typename F>
void WBcoeff<T,F>::extend_w(std::size_t m, std::size_t N) {

  std::size_t Nw, Nj, M = N+1;
  
  if (m+1 > w.size())
    w.resize(m+1);
  
  Nw = w[m].size();

  std::size_t i;

  if (N > Nw) {
    Nj = (N-Nw)/(Nw+1)+1;
    w[m].resize(N);
    i = N;
    F weights;
    if (Nw) {
      for (auto rwi = w[m].rbegin(); rwi != w[m].rend(); ++rwi, --i) {
	if (i%Nj) {
	  *rwi = weights(T(i)/T(M)); //exp(-T(M)*T(M)/(T(i)*T(M-i)));
	} else {
	  *rwi = w[m][i/Nj-1];
	}
      }
    } else {
      for (auto rwi = w[m].rbegin(); rwi != w[m].rend(); ++rwi, --i) {
	*rwi = weights(T(i)/T(M)); //exp(-T(M)*T(M)/T(i)*T(M-i));
	
      }
    }
  }
}


template<typename T, typename F, typename IT>
T WBsum(WBcoeff<T,F>& coeff, IT ui, std::size_t N) {
  
  std::size_t m, n, Nj;
  coeff.find(m,n,N);
  coeff.extend_w(m,N);
  Nj = (coeff.sizew(m)-N)/(N+1)+1;
  
  IT rui = ui + (N-1);
  T sum1(0), sum2(0);
  typename std::vector<T>::const_iterator wi = coeff.cbeginw(m) + (Nj-1);

  coeff.alloc_A(m,n);
  if (coeff.getA(m,n) == T(0)) {
    T A(0);
    for (std::size_t i = 0; i!=N/2; ++i, ++ui, --rui, wi+=Nj) {
      A += *wi;
      sum1 += (*wi)*(*ui);
      sum2 += (*wi)*(*rui);
    }
    A *= 2;
    if ( N%2 == 1) {
      A += *wi;
      sum1 += (*wi)*(*ui);
    }
    sum1 += sum2;
    coeff.setA(m,n,A);
  } else {
    for (std::size_t i = 0; i!=N/2; ++i, ++ui, --rui, wi+=Nj) {
      sum1 += (*wi)*(*ui);
      sum2 += (*wi)*(*rui);
    }
    if ( N%2 == 1 ) {
      sum1 += (*wi)*(*ui);
    }
    sum1 += sum2;
  }

  return sum1/coeff.getA(m,n);
}

template<typename T, typename F, typename IT>
T WBsum_old(WBcoeff<T,F>& coeff, IT ui, std::size_t N) {
 
  std::size_t m, n, Nj;
  coeff.find(m,n,N);
  coeff.extend_w(m,N);
  Nj = (coeff.sizew(m)-N)/(N+1)+1;
  
  T sum(0);
  typename std::vector<T>::const_iterator wi = coeff.cbeginw(m) + (Nj-1);

  coeff.alloc_A(m,n);
  if (coeff.getA(m,n) == T(0)) {
    T A(0);
    for (std::size_t i = 0; i!=N; ++i, ++ui, wi+=Nj) {
      A += *wi;
      sum += (*wi)*(*ui);
    }
    coeff.setA(m,n,A);
  } else {
    for (std::size_t i = 0; i!=N; ++i, ++ui, wi+=Nj) {
      sum += (*wi)*(*ui);
    }
  }
  
  return sum/coeff.getA(m,n);
}

template<typename T, typename F, typename IT>
T WBKahansum(WBcoeff<T,F>& coeff, IT ui, std::size_t N) {
 
  std::size_t m, n, Nj;
  coeff.find(m,n,N);
  coeff.extend_w(m,N);
  Nj = (coeff.sizew(m)-N)/(N+1)+1;
  
  T sum(0), csum(0), cA(0), y, t;
  typename std::vector<T>::const_iterator wi = coeff.cbeginw(m) + (Nj-1);

  coeff.alloc_A(m,n);
  if (coeff.getA(m,n) == T(0)) {
    T A(0);
    for (std::size_t i = 0; i!=N; ++i, ++ui, wi+=Nj) {
      y = *wi - cA;
      t = A + y;
      cA = (t - A) - y;
      A = t;
      y = (*wi)*(*ui) - csum;
      t = sum + y;
      csum = (t - sum) - y;
      sum = t;
    }
    coeff.setA(m,n,A);
  } else {
    for (std::size_t i = 0; i!=N; ++i, ++ui, wi+=Nj) {
      y = (*wi)*(*ui) - csum;
      t = sum + y;
      csum = (t - sum) - y;
      sum = t;
    }
  }
  
  return sum/coeff.getA(m,n);
}

#endif
