
#ifndef _TRANSFORM_TOOLS_
#define _TRANSFORM_TOOLS_
//#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

#include <cmath>
#include <vector>
#include <iterator>
#include "LinAlgTools.hpp"


template<typename Iterator>
void _FOURIER_TRANSFORM_POWER_OF_TWO_(Iterator data, std::size_t nn, int isign);
template<typename ValueType, typename Iterator, typename CIterator>
void _FQFT_COMPLEX_TWO_SIDED_(Iterator output, CIterator input, std::size_t n, ValueType omega);
template<typename ScalarType, typename Iterator, typename CIterator>
void _FQFT_COMPLEX_TWO_SIDED_VEC_(Iterator output, CIterator input, std::size_t n, ScalarType omega);
template<typename ScalarType, typename Iterator, typename CIterator>
void _FQFT_COMPLEX_TWO_SIDED_VEC_V2_(Iterator output, const std::size_t m, CIterator input, const std::size_t n, ScalarType omega);

template<typename Array>
inline void FFT(Array& v);
template<typename Array>
inline void IFFT(Array& v);

template<typename ValueType, typename Array>
inline void QFT(Array& v, ValueType omega);

template<typename Array>
inline void fourier_coeff_ug(Array& v);
template<typename Array>
inline void fourier_series_ug(Array& v);

template<typename ValueType, typename Array>
inline void fourier_coeff_qg(Array& out, const Array& in, ValueType omega);
template<typename ValueType, typename Array>
inline void fourier_coeff_real_qg(Array& out, const Array& in, ValueType omega);

template<typename ValueType, typename Array>
void phase_translate_fourier(Array& data, ValueType w);

template<typename ScalarType, typename Array>
ScalarType fourier_1_norm(const Array& data, const ScalarType& rho);
template<typename ScalarType, typename Array>
ScalarType fourier_max_norm(const Array& data, const ScalarType& rho);
template<typename ScalarType, typename Array>
ScalarType fourier_padding_max_norm(const Array& data, const ScalarType& rho);

template<typename ScalarType>
inline ScalarType sqr(ScalarType x)
{
    return x*x;
}

template<typename ScalarType>
ScalarType _FOURIER_NORM_1_(const ScalarType& _REAL_, const ScalarType& _IMAG_)
{
    using std::sqrt;
    return sqrt(sqr(_REAL_)+sqr(_IMAG_));
}

template<typename ScalarType, std::size_t dim>
ScalarType _FOURIER_NORM_1_(const Vec<ScalarType,dim>& _REAL_, const Vec<ScalarType,dim>& _IMAG_)
{
    using std::sqrt;
    ScalarType res(0);
    for (std::size_t i=0; i!=dim; ++i)
        res += sqrt(sqr(_REAL_[i])+sqr(_IMAG_[i]));
    return res;
}

template<typename ScalarType>
ScalarType _FOURIER_NORM_MAX_(const ScalarType& _REAL_, const ScalarType& _IMAG_)
{
    using std::sqrt;
    return sqrt(sqr(_REAL_)+sqr(_IMAG_));
}

template<typename ScalarType, std::size_t dim>
ScalarType _FOURIER_NORM_MAX_(const Vec<ScalarType,dim>& _REAL_, const Vec<ScalarType,dim>& _IMAG_)
{
    using std::sqrt;
    using std::max;
    ScalarType res(0), temp(0);
    for (std::size_t i=0; i!=dim; ++i) {
        temp = sqrt(sqr(_REAL_[i])+sqr(_IMAG_[i]));
        res = max(res,temp);
    }
    return res;
}
template<typename ScalarType, std::size_t dim>
ScalarType _FOURIER_NORM_MAX_TEST_(const Vec<ScalarType,dim>& _REAL_, const Vec<ScalarType,dim>& _IMAG_)
{
    using std::sqrt;
    using std::max;
    ScalarType res(0), temp(0);
    for (std::size_t i=0; i!=dim; ++i) {
        temp = sqrt(sqr(_REAL_[i])+sqr(_IMAG_[i]));
        res = max(res,temp);
        std::cout << "temp = " << temp << "\nres=" << res << std::endl;
    }
    return res;
}

template<typename _vector_>
struct linear_traits {
    using scalar_type = typename _vector_::value_type;
    static constexpr std::size_t dim = _vector_::dimension;
};
template<>
struct linear_traits<double> {
    using scalar_type = double;
    static constexpr std::size_t dim = 1;
};
template<>
struct linear_traits<my_mpfr> {
    using scalar_type = my_mpfr;
    static constexpr std::size_t dim = 1;
};
template<>
struct linear_traits<INTERVAL> {
    using scalar_type = INTERVAL;
    static constexpr std::size_t dim = 1;
};


template<typename Array>
inline void FFT(Array& v) {
    _FOURIER_TRANSFORM_POWER_OF_TWO_(v.begin(), v.size()/2, 1);
}

template<typename Array>
inline void IFFT(Array& v) {
    //using ValueType = typename Array::value_type;
    using ScalarType = typename linear_traits<typename Array::value_type>::scalar_type;
    _FOURIER_TRANSFORM_POWER_OF_TWO_(v.begin(), v.size()/2, -1);
    ScalarType n = v.size()/2;
    for(auto& vi: v)
        vi /= n;
}

//template<typename ValueType, typename Array>
//inline void QFT(Array& f, ValueType omega) {
//    _FQFT_COMPLEX_TWO_SIDED_VEC_(f.begin(), f.size()/2, omega);
//}

template<typename Array>
inline void fourier_coeff_ug(Array& f) {
    IFFT(f);
}
template<typename Array>
inline void fourier_series_ug(Array& f) {
    FFT(f);
}

template<typename ValueType, typename Array>
inline void QFT(Array& out, const Array& in, ValueType omega) {
    //out.resize(in.size()-2);
    //_FQFT_COMPLEX_TWO_SIDED_VEC_(out.begin(), in.cbegin(), in.size()/2, omega);
    _FQFT_COMPLEX_TWO_SIDED_VEC_V2_(out.begin(), out.size()/2, in.cbegin(), in.size()/2, omega);
}

template<typename ScalarType, typename Array>
inline void fourier_coeff_qg(Array& out, const Array& in, ScalarType omega) {
    auto weights = [](ScalarType t) {
        using std::exp;
        return exp(ScalarType(1.0)/(t*(t-ScalarType(1.0))));
    };
    using ValueType = typename Array::value_type;
    std::vector<ValueType> win(in.size(),ValueType(0));
    using std::exp;
    auto ti = win.begin();
    auto ii = in.cbegin();
    std::size_t n = in.size()/2;
    ScalarType A(0.0), temp;
    for (std::size_t j=0; j!=n; ++j) {
        temp = weights(ScalarType(j+1)/ScalarType(n+1));
        A += temp;
        *ti = *ii*temp;
        ++ti; ++ii;
        *ti = *ii*temp;
        ++ti; ++ii;
    }

    QFT(out,win,-omega);

    for (auto& o: out)
        o /= A;
}

//template<typename ValueType, typename Array>
//inline void fourier_coeff_real_qg(Array& out, const Array& in, ValueType omega) {
//    std::vector<ValueType> temp(2*in.size(),ValueType(0));
//    auto ti = temp.begin();
//    auto ii = in.cbegin();
//    while (ii != in.end()) {
//        *ti = *ii;
//        ti += 2; ++ii;
//    }
//    QFT(out,temp,-omega);
//}

template<typename Array>
INTERVAL fourier_series_real_range(const Array& f, std::size_t N) {
    //using ValueType = typename Array::value_type;
    const INTERVAL twopi = INTERVAL(2)*pi<INTERVAL>();
    std::vector<INTERVAL> f_ext(2*N,INTERVAL(0));
    std::size_t Nf = f.size()/4;
    f_ext[0] = f[0];
    f_ext[1] = f[1];
    for (std::size_t k=1; k!=Nf; ++k) {
        INTERVAL cosval = cos(twopi*k*INTERVAL(0,1)/N);
        INTERVAL sinval = sin(twopi*k*INTERVAL(0,1)/N);
        //std::cout << "cosval = " << cosval << "\nsinval = " << sinval << "\n";
        f_ext[2*k] = f[2*k]*cosval - f[2*k+1]*sinval;
        f_ext[2*k+1] = f[2*k+1]*cosval + f[2*k]*sinval;
        f_ext[2*N-2*k] = f[4*Nf-2*k]*cosval + f[4*Nf-2*k+1]*sinval;
        f_ext[2*N-2*k+1] = f[4*Nf-2*k+1]*cosval - f[4*Nf-2*k]*sinval;
    }
    FFT(f_ext);
    INTERVAL frange(f_ext[0]);
    for (std::size_t j=1; j!=N; ++j) {
        frange = convhull(frange,f_ext[2*j]);
    }

    return frange;
}

template<typename Array>
INTERVAL fourier_series_maxnorm_range_test(const Array& f, INTERVAL rho, std::size_t N) {
    using ValueType = typename Array::value_type;
    const INTERVAL twopi = INTERVAL(2)*pi<INTERVAL>();
    std::vector<ValueType> f_ext(2*N,ValueType(0));
    std::size_t Nf = f.size()/4;
    f_ext[0] = f[0];
    f_ext[1] = f[1];
    //std::cout << "Skip using rho = " << rho << "\n";
    std::size_t nonzero = 0;
    for (std::size_t k=1; k!=Nf; ++k) {
        INTERVAL expval = exp(twopi*k*INTERVAL(rho.right_l()));
        //std::cout << "expval = " << expval << "\n";
        INTERVAL cosval = cos(twopi*k*INTERVAL(0,1)/N);
        INTERVAL sinval = sin(twopi*k*INTERVAL(0,1)/N);
        //std::cout << "cosval = " << cosval << "\nsinval = " << sinval << "\n";
        f_ext[2*k] = expval*(f[2*k]*cosval - f[2*k+1]*sinval);
        f_ext[2*k+1] = expval*(f[2*k+1]*cosval + f[2*k]*sinval);
        f_ext[N-2*k] = expval*(f[2*Nf-2*k]*cosval + f[2*Nf-2*k+1]*sinval);
        f_ext[N-2*k+1] = expval*(f[2*Nf-2*k+1]*cosval - f[2*Nf-2*k]*sinval);
        if( _FOURIER_NORM_MAX_(f_ext[2*k],f_ext[2*k+1]).right_l() != my_mpfr(0) ) {
            ++nonzero;
            //std::cout << "expval = " << expval << "\n";
            //std::cout << "|f| = " << _FOURIER_NORM_MAX_(f[2*k],f[2*k+1]) << "\n";
            //std::cout << "|f_ext| = " << _FOURIER_NORM_MAX_(f_ext[2*k],f_ext[2*k+1]) << "\n";
        }
    }
    std::cout << "Nf = " << Nf << ", nonzero = " << nonzero << std::endl;
    FFT(f_ext);
    INTERVAL frange(_FOURIER_NORM_MAX_(f_ext[0],f_ext[1]));
    std::cout << "|f_ext[0]| = " << _FOURIER_NORM_MAX_(f_ext[0],f_ext[1]) << "\n";
    for (std::size_t j=1; j!=N; ++j) {
        frange = convhull(frange,_FOURIER_NORM_MAX_(f_ext[2*j],f_ext[2*j+1]));
    }
    return frange;
}

template<typename Array>
INTERVAL fourier_series_maxnorm_range(const Array& f, INTERVAL rho, std::size_t N) {
    using ValueType = typename Array::value_type;
    const INTERVAL twopi = INTERVAL(2)*pi<INTERVAL>();
    std::vector<ValueType> f_ext(2*N,ValueType(0));
    std::size_t Nf = f.size()/4;
    f_ext[0] = f[0];
    f_ext[1] = f[1];
    for (std::size_t k=1; k!=Nf; ++k) {
        INTERVAL cosval = cos(twopi*k*INTERVAL(0,1)/N);
        INTERVAL sinval = sin(twopi*k*INTERVAL(0,1)/N);
        //std::cout << "cosval = " << cosval << "\nsinval = " << sinval << "\n";
        f_ext[2*k] = f[2*k]*cosval - f[2*k+1]*sinval;
        f_ext[2*k+1] = f[2*k+1]*cosval + f[2*k]*sinval;
        f_ext[N-2*k] = f[2*Nf-2*k]*cosval + f[2*Nf-2*k+1]*sinval;
        f_ext[N-2*k+1] = f[2*Nf-2*k+1]*cosval - f[2*Nf-2*k]*sinval;
    }
    FFT(f_ext);
    INTERVAL frange(_FOURIER_NORM_MAX_(f_ext[0],f_ext[1]));
    //std::cout << "|f_ext[0]| = " << _FOURIER_NORM_MAX_(f_ext[0],f_ext[1]) << "\n";
    for (std::size_t j=1; j!=N; ++j) {
        frange = convhull(frange,_FOURIER_NORM_MAX_(f_ext[2*j],f_ext[2*j+1]));
    }
    if (rho.left_l() != my_mpfr(0) || rho.right_l() != my_mpfr(0)) {
        //INTERVAL padding = fourier_max_norm(f,rho) - fourier_max_norm(f,INTERVAL(0));
        //std::cout << "padding = " << padding << std::endl;
        INTERVAL padding = fourier_padding_max_norm(f,rho);
        //std::cout << "padding = " << padding << std::endl;
        padding *= INTERVAL(-1.0,1.0);
        frange += padding;
    }
    return frange;
}

template<typename ScalarType, typename Array>
void phase_translate_fourier(Array& data, ScalarType phase) {
    using std::sin;
    using std::cos;
    using ValueType = typename Array::value_type;

    std::size_t n = data.size()/2;

    ScalarType pi2 = ScalarType(2.0)*pi<ScalarType>();
    ScalarType tempr, tempi;
    ValueType datar;

    for (std::size_t k=1; k!=n/2; k++) {
        //std::cout << "\n k = " << k << std::endl;

        //std::cout << "(data["<<2*k<<"],data["<<2*k+1<<"]) = ("<<data[2*k]<<","<<data[2*k+1]<<")\n";
        //std::cout << "(data["<<2*(n-k)<<"],data["<<2*(n-k)+1<<"]) = ("<<data[2*(n-k)]<<","<<data[2*(n-k)+1]<<")\n";

        tempr = cos(pi2*phase*k);
        tempi = sin(pi2*phase*k);

        //std::cout << "tempr = " << tempr << ", tempi = " << tempi << "\n";

        //std::cout << "data["<<2*k<<"] = "<<data[2*k]<<"*"<<tempr<<" - "<<data[2*k+1]<<"*"<<tempi<<"\n";
        //std::cout << "data["<<2*k+1<<"] = "<<data[2*k]<<"*"<<tempi<<" + "<<data[2*k+1]<<"*"<<tempr<<"\n";

        //std::cout << "data["<<2*(n-k)<<"] = "<<data[2*(n-k)]<<"*"<<tempr<<" + "<<data[2*(n-k)+1]<<"*"<<tempi<<"\n";
        //std::cout << "data["<<2*(n-k)+1<<"] = -"<<data[2*(n-k)]<<"*"<<tempi<<" + "<<data[2*(n-k)+1]<<"*"<<tempr<<"\n";

        datar = data[2*k];
        data[2*k] = datar*tempr - data[2*k+1]*tempi;
        data[2*k+1] = datar*tempi + data[2*k+1]*tempr;

        datar = data[2*(n-k)];
        data[2*(n-k)] = datar*tempr + data[2*(n-k)+1]*tempi;
        data[2*(n-k)+1] = -datar*tempi + data[2*(n-k)+1]*tempr;

        //std::cout << "(data["<<2*k<<"],data["<<2*k+1<<"]) = ("<<data[2*k]<<","<<data[2*k+1]<<")\n";
        //std::cout << "(data["<<2*(n-k)<<"],data["<<2*(n-k)+1<<"]) = ("<<data[2*(n-k)]<<","<<data[2*(n-k)+1]<<")\n";
    }

    //std::cout << "\n k = " << n/2 << std::endl;

    //std::cout << "(data["<<n/2<<"],data["<<n/2+1<<"]) = ("<<data[n]<<","<<data[n+1]<<")\n";

    tempr = cos(pi2*phase*n/2);
    tempi = sin(pi2*phase*n/2);

    //std::cout << "tempr = " << tempr << ", tempi = " << tempi << "\n";

    //std::cout << "data["<<n<<"] = "<<data[n]<<"*"<<tempr<<" - "<<data[n+1]<<"*"<<tempi<<"\n";
    //std::cout << "data["<<n+1<<"] = "<<data[n]<<"*"<<tempi<<" + "<<data[n+1]<<"*"<<tempr<<"\n";

    datar = data[n];
    data[n] = datar*tempr + data[n+1]*tempi;
    data[n+1] = -datar*tempi + data[n+1]*tempr;

    //std::cout << "(data["<<n/2<<"],data["<<n/2+1<<"]) = ("<<data[n]<<","<<data[n+1]<<")\n" << std::endl;

}

template<typename ScalarType, typename Array>
ScalarType fourier_1_norm(const Array& data, const ScalarType& rho)
{
    //using ValueType = typename Array::value_type;
    using std::exp;
    std::size_t n = data.size()/2, n2 = n/2;
    ScalarType res(0), weight(0), pi2 = 2*pi<ScalarType>();
    for (std::size_t k=0; k!= n2; ++k) {
        weight = exp(pi2*rho*k);
        res += _FOURIER_NORM_1_(data[2*k],data[2*k+1])*weight;
    }
    for (std::size_t k=1; k!= n2+1; ++k) {
        weight = exp(pi2*rho*k);
        res += _FOURIER_NORM_1_(data[2*(n-k)],data[2*(n-k)+1])*weight;
    }
    return res;
}

template<typename ScalarType, typename Array>
ScalarType fourier_max_norm(const Array& data, const ScalarType& rho)
{
    //using ValueType = typename Array::value_type;
    using std::exp;
    std::size_t n = data.size()/2, n2 = n/2;
    ScalarType res(0), weight(0), pi2 = 2*pi<ScalarType>();
    for (std::size_t k=0; k!= n2; ++k) {
        weight = exp(pi2*rho*k);
        res += _FOURIER_NORM_MAX_(data[2*k],data[2*k+1])*weight;
    }
    for (std::size_t k=1; k!= n2+1; ++k) {
        weight = exp(pi2*rho*k);
        res += _FOURIER_NORM_MAX_(data[2*(n-k)],data[2*(n-k)+1])*weight;
    }
    return res;
}

template<typename ScalarType, typename Array>
ScalarType fourier_padding_max_norm(const Array& data, const ScalarType& rho)
{
    //using ValueType = typename Array::value_type;
    using std::exp;
    using std::abs;
    std::size_t n = data.size()/2, n2 = n/2;
    ScalarType res1(0), res2(0), normfront(0), normback(0), weight(0), pi2 = 2*pi<ScalarType>();
    for (std::size_t k=1; k!= n2; ++k) {
        normfront = _FOURIER_NORM_MAX_(data[2*k],data[2*k+1]);
        normback = _FOURIER_NORM_MAX_(data[2*(n-k)],data[2*(n-k)+1]);
        weight = abs(exp(pi2*rho*k)-ScalarType(1));
        res1 += normfront*weight;
        res2 += normback*weight;
        weight = abs(exp(-pi2*rho*k)-ScalarType(1));
        res1 += normback*weight;
        res2 += normfront*weight;
    }
    normfront = _FOURIER_NORM_MAX_(data[2*n2],data[2*n2+1]);
    weight = exp(-pi2*rho*n2)-ScalarType(1);
    res1 += normfront*weight;
    weight = exp(pi2*rho*n2)-ScalarType(1);
    res2 += normfront*weight;
    return max(res1,res2);
}

template<typename ScalarType>
ScalarType DFT_error_const(const ScalarType rho, ScalarType rho_hat, std::size_t NF)
{
    using std::exp;
    ScalarType S1, S2, T, exponential;
    const ScalarType pi1 = pi<ScalarType>(), pi2 = 2*pi1;

    exponential = exp(-pi2*rho_hat*NF);
    S1 = S2 = exponential/(1 - exponential);

    exponential = exp(-pi2*(rho_hat+rho));
    S1 *= (exponential + 1)/(exponential - 1)*(1 - exp(pi1*(rho_hat+rho)*NF));

    exponential = exp(pi2*(rho_hat-rho));
    T = (exponential + 1)/(exponential - 1);
    S2 *= T;

    exponential = exp(-pi1*(rho_hat-rho)*NF);
    S2 *= 1 - exponential;
    T *= exponential;

    return S1 + S2 + T;
}

template<typename ScalarType>
ScalarType DFT_error_coeff(std::size_t k, ScalarType rho, std::size_t NF)
{
    using std::exp;
    ScalarType S, exponential;
    const ScalarType pi2 = 2*pi<ScalarType>();

    exponential = exp(-pi2*rho*NF);
    S = exponential/(1 - exponential);

    if (k == 0)
        exponential = 2;
    else
        exponential = exp(pi2*rho*k) + exp(-pi2*rho*k);

    return S*exponential;
}

template<typename Iterator>
void _FOURIER_TRANSFORM_POWER_OF_TWO_(Iterator data, std::size_t nn, int isign)
{
    using ValueType = typename std::iterator_traits<Iterator>::value_type;
    using ScalarType = typename linear_traits<ValueType>::scalar_type;

    //std::cout << "Starting FFT" << std::endl;
    --data; // Decrement pointer for zero-offset data array

    std::size_t n,mmax,m,j,istep,i;
    ScalarType wtemp,wr,wpr,wpi,wi,theta;
    ValueType tempr,tempi;
    ScalarType twoPI = 2*pi<ScalarType>();
    n=nn << 1; // nn - complex length, n - real length
    j=1;

    //std::cout << "nn = " << nn << ", n = " << n << std::endl;
    //std::cout << "data = \n";
    //for (std::size_t i=1; i!=n+1; i+=2)
    //    std::cout << "(" << data[i] << "," << data[i+1] << ")\n";

    //std::cout << "Do bit reversal\n" << std::endl;
    for (i=1;i<n;i+=2) {
        //std::cout << "i = " << i << " j = " << j << std::endl;
        if (j > i) {
            using std::swap;
            swap(data[j],data[i]);
            swap(data[j+1],data[i+1]);
        }
        m=nn;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    //std::cout << "data =  \n";
    //for (std::size_t i=1; i!=n+1; i+=2)
    //    std::cout << "(" << data[i] << "," << data[i+1] << ")\n";


    //std::cout << "Do transformation" << std::endl;
    mmax=2;
    while (n > mmax) {
        istep = mmax << 1;
        theta = isign*(twoPI/mmax);
        wtemp = sin(ScalarType(0.5)*theta);
        wpr = ScalarType(-2.0)*wtemp*wtemp;
        wpi = sin(theta);
        wr = ScalarType(1.0);
        wi = ScalarType(0.0);
        //std::cout << "Start loop" << std::endl;
        for (std::size_t m=1;m<mmax;m+=2) {
            //std::cout << "m=" << m << std::endl;
            for (std::size_t i=m;i<=n;i+=istep) {
                //std::cout << "i=" << i << std::endl;
                j=i+mmax;
                //std::cout << "j=" << j << std::endl;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i]=data[i]+tempr;
                data[i+1]=data[i+1]+tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;

        //std::cout << "data =  \n";
        //for (std::size_t i=1; i!=n+1; i+=2)
        //    std::cout << "(" << data[i] << "," << data[i+1] << ")\n";

    }
    //std::cout << "FFT finished\n" << std::endl;
}

template<typename ValueType>
void _FOURIER_TRANSFORM_TWO_REAL_FUNCTIONS_POWER_OF_TWO_(ValueType data1[], ValueType data2[], ValueType fft1[], ValueType fft2[], std::size_t n) {

    --data1; --data2; --fft1; --fft2; // decrement pointers

    std::size_t nn3,nn2,jj,j;
    ValueType rep,rem,aip,aim;
    nn3=1+(nn2=2+n+n);
    for (j=1,jj=2;j<=n;j++,jj+=2) {
        fft1[jj-1]=data1[j];
        fft1[jj]=data2[j];
    }
    _FOURIER_TRANSFORM_POWER_OF_TWO_(fft1+1,n,1); // Remeber to undecrement pointer
    fft2[1]=fft1[2];
    fft1[2]=fft2[2]=0.0;
    for (j=3;j<=n+1;j+=2) {
        rep=0.5*(fft1[j]+fft1[nn2-j]);
        rem=0.5*(fft1[j]-fft1[nn2-j]);
        aip=0.5*(fft1[j+1]+fft1[nn3-j]);
        aim=0.5*(fft1[j+1]-fft1[nn3-j]);
        fft1[j]=rep;
        fft1[j+1]=aim;
        fft1[nn2-j]=rep;
        fft1[nn3-j] = -aim;
        fft2[j]=aip;
        fft2[j+1] = -rem;
        fft2[nn2-j]=aip;
        fft2[nn3-j]=rem;
    }
}


// Note that n here is complex length of input
// Consider removing the calls to cos and sin and calculate the values recursively instead
template<typename ValueType, typename Iterator, typename CIterator>
void _FQFT_COMPLEX_TWO_SIDED_(Iterator output, CIterator input, std::size_t n, ValueType omega) {

    static_assert(std::is_same< ValueType, typename std::iterator_traits<Iterator>::value_type >::value,"Iterator:value_type must be ValueType");
    static_assert(std::is_same< ValueType, typename std::iterator_traits<CIterator>::value_type >::value,"CIterator:value_type must be ValueType");

    using std::sin; // For machine types
    using std::cos;
    ValueType PI = pi<ValueType>();

    //ValueType* print_ptr = input;
    //for (std::size_t j=0; j!=n; ++j, print_ptr+=2)
    //    std::cout << "input["<<j<<"] = ("<<*print_ptr<<","<<*(print_ptr+1)<<")\n";
    //print_ptr = output;
    //for (std::size_t j=0; j!=n-1; ++j, print_ptr+=2)
    //    std::cout << "output["<<j<<"] = ("<<*print_ptr<<","<<*(print_ptr+1)<<")\n";
    //std::cout << "omega = " << omega << std::endl;

    // For debugging
    //using std::cout; using std::endl; //complex<ValueType> temp_print(0.0);

    // n = size of input
    const std::size_t m = n-1, m2 = m/2, N = 2*m, N2 = N/2; // NOTE: n = 2*m2 + 1! since n is odd
    //output.resize(N);
    std::vector<ValueType> g(2*N,ValueType(0.0)); 
    auto gptf = g.begin();
    auto gptb = g.end()-2; // Points to real part of last complex element
    std::vector<ValueType> h(2*N,ValueType(0.0)); 
    auto hptf = h.begin();
    auto hptb = h.end()-2; // Points to real part of last complex element

    //std::cout << "Starting _FQFT_POWER_OF_TWO_PLUS_ONE_" << std::endl;
    //std::cout <<" with n = " << n << ", m2 = " << m2 << ", N = " << N << ", N2 = " << N2 << std::endl;

    //-----------Construct inputs for convolution-------------------
    CIterator iptf = input+(2*m2);
    CIterator iptb = iptf-2;
    Iterator optf = output;
    Iterator optb = output+2*(m-1);
    
    *hptf = ValueType(1.0);
    ++hptf;
    *hptf = ValueType(0.0);
    ++hptf;

    optf += 2; // Skip saving this weigth
    
    *gptf = *iptf;
    ++gptf; ++iptf;
    *gptf = *iptf;
    ++gptf; ++iptf;

    //cout << "h[" << 0 << "] = (" << *(hptf-2) << "," << *(hptf-1) << ")\n";
    //cout << "g[" << 0 << "] = (" << *(gptf-2) << "," << *(gptf-1) << ")\n";

    for (std::size_t j=1; j!=m2+1; ++j) {
        *hptf =  cos(PI*(j*j*omega));
        *hptb = *hptf;
        *(hptf+1) = -sin(PI*(j*j*omega));
        *(hptb+1) = *(hptf+1);

        *optf = *hptf; // Save weights for output
        ++optf; 
        *optf = -*(hptf+1);
        ++optf;

        *gptf = (*hptf)*(*iptf) + (*(hptf+1))*(*(iptf+1));
        ++gptf;
        *gptf = (*hptf)*(*(iptf+1)) - (*(hptf+1))*(*iptf);
        ++gptf;
        hptf += 2; iptf += 2;
        
        *gptb = (*hptb)*(*iptb) + (*(hptb+1))*(*(iptb+1));
        *(gptb+1) = (*hptb)*(*(iptb+1)) - (*(hptb+1))*(*iptb);
        gptb -= 2; hptb -= 2; iptb -= 2;


        //cout << "h["<<j<<"] = (" << *(hptf-2) << "," << *(hptf-1) << "), h["<<N-j<<"] = (" << *(hptb+2) << "," << *(hptb+3) << ")\n";
        //cout << "g["<<j<<"] = (" << *(gptf-2) << "," << *(gptf-1) << "), g["<<N-j<<"] = (" << *(gptb+2) << "," << *(gptb+3) << ")\n";
    }
    for (std::size_t j=m2+1; j!=N2; ++j) {
        *hptf =  cos(PI*(j*j*omega));
        *hptb = *hptf;
        *(hptf+1) = -sin(PI*(j*j*omega));
        *(hptb+1) = *(hptf+1);
        hptf += 2; hptb -= 2;

        gptf += 2; gptb -= 2;
        //cout << "h["<<j<<"] = (" << *(hptf-2) << "," << *(hptf-1) << "), h["<<N-j<<"] = (" << *(hptb+2) << "," << *(hptb+3) << ")\n";
        //cout << "g["<<j<<"] = (" << *(gptf-2) << "," << *(gptf-1) << "), g["<<N-j<<"] = (" << *(gptb+2) << "," << *(gptb+3) << ")\n";
    }

    *hptf = cos(PI*(N2*N2*omega));
    *(hptf+1) = -sin(PI*(N2*N2*omega));

    //cout << "h["<<N2<<"] = (" << *hptf << "," << *(hptf+1) << ")\n";
    //cout << "g["<<N2<<"] = (" << *gptf << "," << *(gptf+1) << ")\n";

    //for (std::size_t j=0; j!=N; ++j) 
    //    cout << "h["<<j<<"] = ("<< h[2*j] <<","<< h[2*j+1] <<")\n";
    //for (std::size_t j=0; j!=N; ++j) 
    //    cout << "g["<<j<<"] = ("<< g[2*j] <<","<< g[2*j+1] <<")\n";
    

    //Do convolution using Fourier transform

    _FOURIER_TRANSFORM_POWER_OF_TWO_((ValueType*)g.data(), N, 1);
    _FOURIER_TRANSFORM_POWER_OF_TWO_((ValueType*)h.data(), N, 1);

    for (std::size_t j=0; j!=N; ++j) {
        ValueType mtemp = h[2*j];
        h[2*j] = g[2*j]*mtemp - g[2*j+1]*h[2*j+1];
        h[2*j+1] = g[2*j]*h[2*j+1] + g[2*j+1]*mtemp;
    }

    _FOURIER_TRANSFORM_POWER_OF_TWO_((ValueType*)h.data(), N, -1);

    for(auto& vi: h)
        vi /= (ValueType) N; 

    //for (std::size_t j=0; j!=N; ++j) 
    //    cout << "h["<<j<<"] = ("<< h[2*j] <<","<< h[2*j+1] <<")\n";

    //--------------------------------------------------------------------------------------//

    // ----- Resize output from N to m = N/2 and multiply by weight factors ------
    ValueType rtemp, itemp;// , mtemp;
    optf = output;
    hptf = h.begin();
    hptb = h.end()-2;

    // First element has weight 1
    *optf = *hptf; // Real part
    ++optf; ++hptf;
    *optf = *hptf; // Imaginary part
    ++optf; ++hptf;
    
    for (std::size_t k=1; k!=m2; ++k) {
        rtemp = *optf;     //cos(PI*(k*k*omega));
        itemp = *(optf+1); //sin(PI*(k*k*omega));

        *optf = rtemp*(*hptf) - itemp*(*(hptf+1));
        *(optf+1) = rtemp*(*(hptf+1)) + itemp*(*hptf);
        optf += 2; hptf += 2;

        *optb = rtemp*(*hptb) - itemp*(*(hptb+1));
        *(optb+1) = rtemp*(*(hptb+1)) + itemp*(*hptb);
        optb -= 2; hptb -= 2;
    }

    rtemp = cos(PI*(m2*m2*omega));
    itemp = sin(PI*(m2*m2*omega));
    *optb = rtemp*(*hptb) - itemp*(*(hptb+1));
    *(optb+1) = rtemp*(*(hptb+1)) + itemp*(*hptb);

    //print_ptr = output;
    //for (std::size_t j=0; j!=n-1; ++j, print_ptr+=2)
    //    std::cout << "output["<<j<<"] = ("<<*print_ptr<<","<<*(print_ptr+1)<<")\n";

    //std::cout << "_FQFT_POWER_OF_TWO_PLUS_ONE_ finished" << std::endl;
}

template<typename ScalarType, typename Iterator, typename CIterator>
void _FQFT_COMPLEX_TWO_SIDED_VEC_(Iterator output, CIterator input, std::size_t n, ScalarType omega) {

    static_assert(std::is_same< typename std::iterator_traits<Iterator>::value_type, typename std::iterator_traits<CIterator>::value_type >::value,"Iterator:value_type must be CIterator::value_type");

    using ValueType = typename Iterator::value_type;

    using std::sin; // For machine types
    using std::cos;
    ScalarType PI = pi<ScalarType>();

    //ValueType* print_ptr = input;
    //for (std::size_t j=0; j!=n; ++j, print_ptr+=2)
    //    std::cout << "input["<<j<<"] = ("<<*print_ptr<<","<<*(print_ptr+1)<<")\n";
    //print_ptr = output;
    //for (std::size_t j=0; j!=n-1; ++j, print_ptr+=2)
    //    std::cout << "output["<<j<<"] = ("<<*print_ptr<<","<<*(print_ptr+1)<<")\n";
    //std::cout << "omega = " << omega << std::endl;

    // For debugging
    //using std::cout; using std::endl; //complex<ValueType> temp_print(0.0);

    // n = size of input
    const std::size_t m = n-1, m2 = m/2, N = 2*m, N2 = N/2; // NOTE: n = 2*m2 + 1! since n is odd
    //output.resize(N);
    std::vector<ValueType> g(2*N,ValueType(0.0)); 
    auto gptf = g.begin();
    auto gptb = g.end()-2;
    std::vector<ScalarType> h(2*N,ScalarType(0.0)); 
    auto hptf = h.begin();
    auto hptb = h.end()-2;

    //std::cout << "Starting _FQFT_POWER_OF_TWO_PLUS_ONE_" << std::endl;
    //std::cout <<" with n = " << n << ", m2 = " << m2 << ", N = " << N << ", N2 = " << N2 << std::endl;

    //-----------Construct inputs for convolution-------------------
    CIterator iptf = input+(2*m2);
    CIterator iptb = iptf-2;
    Iterator optf = output;
    Iterator optb = output+2*(m-1);
    
    *hptf = ScalarType(1.0);
    ++hptf;
    *hptf = ScalarType(0.0);
    ++hptf;

    optf += 2; // Skip saving this weigth
    
    *gptf = *iptf;
    ++gptf; ++iptf;
    *gptf = *iptf;
    ++gptf; ++iptf;

    //cout << "h[" << 0 << "] = (" << *(hptf-2) << "," << *(hptf-1) << ")\n";
    //cout << "g[" << 0 << "] = (" << *(gptf-2) << "," << *(gptf-1) << ")\n";

    for (std::size_t j=1; j!=m2+1; ++j) {
        *hptf =  cos(PI*(j*j*omega));
        *hptb = *hptf;
        *(hptf+1) = -sin(PI*(j*j*omega));
        *(hptb+1) = *(hptf+1);

        first_element(*optf) = *hptf; // Save weights for output
        ++optf; 
        first_element(*optf) = -*(hptf+1);
        ++optf;

        *gptf = (*hptf)*(*iptf) + (*(hptf+1))*(*(iptf+1));
        ++gptf;
        *gptf = (*hptf)*(*(iptf+1)) - (*(hptf+1))*(*iptf);
        ++gptf;
        hptf += 2; iptf += 2;
        
        *gptb = (*hptb)*(*iptb) + (*(hptb+1))*(*(iptb+1));
        *(gptb+1) = (*hptb)*(*(iptb+1)) - (*(hptb+1))*(*iptb);
        gptb -= 2; hptb -= 2; iptb -= 2;


        //cout << "h["<<j<<"] = (" << *(hptf-2) << "," << *(hptf-1) << "), h["<<N-j<<"] = (" << *(hptb+2) << "," << *(hptb+3) << ")\n";
        //cout << "g["<<j<<"] = (" << *(gptf-2) << "," << *(gptf-1) << "), g["<<N-j<<"] = (" << *(gptb+2) << "," << *(gptb+3) << ")\n";
    }
    for (std::size_t j=m2+1; j!=N2; ++j) {
        *hptf =  cos(PI*(j*j*omega));
        *hptb = *hptf;
        *(hptf+1) = -sin(PI*(j*j*omega));
        *(hptb+1) = *(hptf+1);
        hptf += 2; hptb -= 2;

        gptf += 2; gptb -= 2;
        //cout << "h["<<j<<"] = (" << *(hptf-2) << "," << *(hptf-1) << "), h["<<N-j<<"] = (" << *(hptb+2) << "," << *(hptb+3) << ")\n";
        //cout << "g["<<j<<"] = (" << *(gptf-2) << "," << *(gptf-1) << "), g["<<N-j<<"] = (" << *(gptb+2) << "," << *(gptb+3) << ")\n";
    }

    *hptf = cos(PI*(N2*N2*omega));
    *(hptf+1) = -sin(PI*(N2*N2*omega));

    //cout << "h["<<N2<<"] = (" << *hptf << "," << *(hptf+1) << ")\n";
    //cout << "g["<<N2<<"] = (" << *gptf << "," << *(gptf+1) << ")\n";

    //for (std::size_t j=0; j!=N; ++j) 
    //    cout << "h["<<j<<"] = ("<< h[2*j] <<","<< h[2*j+1] <<")\n";
    //for (std::size_t j=0; j!=N; ++j) 
    //    cout << "g["<<j<<"] = ("<< g[2*j] <<","<< g[2*j+1] <<")\n";
    

    //Do convolution using Fourier transform

    _FOURIER_TRANSFORM_POWER_OF_TWO_((ValueType*)g.data(), N, 1);
    _FOURIER_TRANSFORM_POWER_OF_TWO_((ScalarType*)h.data(), N, 1);

    for (std::size_t j=0; j!=N; ++j) {
        ValueType mtemp = g[2*j];
        g[2*j] = h[2*j]*mtemp - h[2*j+1]*g[2*j+1];
        g[2*j+1] = h[2*j]*g[2*j+1] + h[2*j+1]*mtemp;
    }

    _FOURIER_TRANSFORM_POWER_OF_TWO_((ValueType*)g.data(), N, -1);

    for(auto& vi: g)
        vi /= (ScalarType) N; 

    //for (std::size_t j=0; j!=N; ++j) 
    //    cout << "h["<<j<<"] = ("<< h[2*j] <<","<< h[2*j+1] <<")\n";

    //--------------------------------------------------------------------------------------//

    // ----- Resize output from N to m = N/2 and multiply by weight factors ------
    ScalarType rtemp, itemp;// , mtemp;
    optf = output;
    gptf = g.begin();
    gptb = g.end()-2;

    // First element has weight 1
    *optf = *gptf; // Real part
    ++optf; ++gptf;
    *optf = *gptf; // Imaginary part
    ++optf; ++gptf;
    
    for (std::size_t k=1; k!=m2; ++k) {
        rtemp = first_element(*optf);     //cos(PI*(k*k*omega));
        itemp = first_element(*(optf+1)); //sin(PI*(k*k*omega));

        *optf = rtemp*(*gptf) - itemp*(*(gptf+1));
        *(optf+1) = rtemp*(*(gptf+1)) + itemp*(*gptf);
        optf += 2; gptf += 2;

        *optb = rtemp*(*gptb) - itemp*(*(gptb+1));
        *(optb+1) = rtemp*(*(gptb+1)) + itemp*(*gptb);
        optb -= 2; gptb -= 2;
    }

    rtemp = cos(PI*(m2*m2*omega));
    itemp = sin(PI*(m2*m2*omega));
    *optb = rtemp*(*gptb) - itemp*(*(gptb+1));
    *(optb+1) = rtemp*(*(gptb+1)) + itemp*(*gptb);

    //print_ptr = output;
    //for (std::size_t j=0; j!=n-1; ++j, print_ptr+=2)
    //    std::cout << "output["<<j<<"] = ("<<*print_ptr<<","<<*(print_ptr+1)<<")\n";

    //std::cout << "_FQFT_POWER_OF_TWO_PLUS_ONE_ finished" << std::endl;
}

// m = complex length of output, n = complex length of input.
// We always assume m <= n
template<typename ScalarType, typename Iterator, typename CIterator>
void _FQFT_COMPLEX_TWO_SIDED_VEC_V2_(Iterator output, const std::size_t m, CIterator input, const std::size_t n, ScalarType omega) {

    static_assert(std::is_same< typename std::iterator_traits<Iterator>::value_type, typename std::iterator_traits<CIterator>::value_type >::value,"Iterator:value_type must be CIterator::value_type");

    using ValueType = typename Iterator::value_type;

    using std::sin; // For machine types
    using std::cos;
    ScalarType PI = pi<ScalarType>();

    //ValueType* print_ptr = input;
    //for (std::size_t j=0; j!=n; ++j, print_ptr+=2)
    //    std::cout << "input["<<j<<"] = ("<<*print_ptr<<","<<*(print_ptr+1)<<")\n";
    //print_ptr = output;
    //for (std::size_t j=0; j!=n-1; ++j, print_ptr+=2)
    //    std::cout << "output["<<j<<"] = ("<<*print_ptr<<","<<*(print_ptr+1)<<")\n";
    //std::cout << "omega = " << omega << std::endl;

    // For debugging
    //using std::cout; using std::endl; //complex<ValueType> temp_print(0.0);

    // m = complex length of output, n = complex length of input.
    const std::size_t m2 = (m+1)/2, n2 = (n+1)/2;
    bool m_is_even = (2*m2 == m), n_is_even = (2*n2 == n);
    std::size_t p = 0, N =1;
    while (n+m > N && p!=53) {
        ++p;
        N = N<<1;
    }
    const std::size_t N2 = N/2;
    //output.resize(N);
    std::vector<ValueType> g(2*N,ValueType(0.0)); 
    auto gptf = g.begin();
    auto gptb = g.end()-2;
    std::vector<ScalarType> h(2*N,ScalarType(0.0)); 
    auto hptf = h.begin();
    auto hptb = h.end()-2;

    //std::cout << "Starting _FQFT_COMPLEX_TWO_SIDED_VEC_V2_" << std::endl;
    //std::cout <<" with n = " << n << ", m = " << m << ", N = " << N << ", N2 = " << N2 << std::endl;

    //-----------Construct inputs for convolution-------------------
    CIterator iptf = input+(2*(n-n2));
    CIterator iptb = iptf-2;
    Iterator optf = output;
    Iterator optb = output+2*(m-1);
    
    *hptf = ScalarType(1.0);
    ++hptf;
    *hptf = ScalarType(0.0);
    ++hptf;

    optf += 2; // Skip saving this weigth
    
    *gptf = *iptf;
    ++gptf; ++iptf;
    *gptf = *iptf;
    ++gptf; ++iptf;

    //cout << "h[" << 0 << "] = (" << *(hptf-2) << "," << *(hptf-1) << ")\n";
    //cout << "g[" << 0 << "] = (" << *(gptf-2) << "," << *(gptf-1) << ")\n";

    for (std::size_t j=1; j!=n2; ++j) {
        //cout << "j = " << j << endl;
        //cout << "g[2*j] = " << g[2*j] << ", g[2*j+1] = " << g[2*j+1] << endl;
        //cout << "h[2*(N-j)] = " << h[2*(N-j)] << ", h[2*(N-j)+1] = " << h[2*(N-j)+1] << endl;

        *hptf =  cos(PI*(j*j*omega));
        *hptb = *hptf;
        *(hptf+1) = -sin(PI*(j*j*omega));
        *(hptb+1) = *(hptf+1);

        if (j<m2) {
            first_element(*optf) = *hptf; // Save weights for output
            ++optf; 
            first_element(*optf) = -*(hptf+1);
            ++optf;
        }
        if ( j == m2 && m_is_even){
            //cout << "j == " << m2 << " && m_is_even" << endl;
            first_element(*optf) = *hptf; // Save weights for output
            ++optf; 
            first_element(*optf) = -*(hptf+1);
            ++optf;
        }


        *gptf = (*hptf)*(*iptf) + (*(hptf+1))*(*(iptf+1));
        ++gptf;
        *gptf = (*hptf)*(*(iptf+1)) - (*(hptf+1))*(*iptf);
        ++gptf;
        hptf += 2; iptf += 2;
        
        *gptb = (*hptb)*(*iptb) + (*(hptb+1))*(*(iptb+1));
        *(gptb+1) = (*hptb)*(*(iptb+1)) - (*(hptb+1))*(*iptb);
        gptb -= 2; hptb -= 2; iptb -= 2;


        //cout << "h["<<j<<"] = (" << *(hptf-2) << "," << *(hptf-1) << "), h["<<N-j<<"] = (" << *(hptb+2) << "," << *(hptb+3) << ")\n";
        //cout << "g["<<j<<"] = (" << *(gptf-2) << "," << *(gptf-1) << "), g["<<N-j<<"] = (" << *(gptb+2) << "," << *(gptb+3) << ")\n";
    }

    //cout << "n2 = " << n2 << endl;

    *hptf =  cos(PI*(n2*n2*omega));
    *hptb = *hptf;
    *(hptf+1) = -sin(PI*(n2*n2*omega));
    *(hptb+1) = *(hptf+1);
    hptf += 2; hptb -= 2;

    if (n_is_even) {
        *gptb = (*hptb)*(*iptb) + (*(hptb+1))*(*(iptb+1));
        *(gptb+1) = (*hptb)*(*(iptb+1)) - (*(hptb+1))*(*iptb);
        gptb -= 2;
    }
    gptf += 2; gptb -= 2;

    //cout << "h["<<n2<<"] = (" << *(hptf-2) << "," << *(hptf-1) << "), h["<<N-n2<<"] = (" << *(hptb+2) << "," << *(hptb+3) << ")\n";
    //cout << "g["<<n2<<"] = (" << *(gptf-2) << "," << *(gptf-1) << "), g["<<N-n2<<"] = (" << *(gptb+2) << "," << *(gptb+3) << ")\n";

    for (std::size_t j=n2+1; j!=N2; ++j) {
        //cout << "j = " << j << endl;

        *hptf =  cos(PI*(j*j*omega));
        *hptb = *hptf;
        *(hptf+1) = -sin(PI*(j*j*omega));
        *(hptb+1) = *(hptf+1);
        hptf += 2; hptb -= 2;

        gptf += 2; gptb -= 2;
        //cout << "h["<<j<<"] = (" << *(hptf-2) << "," << *(hptf-1) << "), h["<<N-j<<"] = (" << *(hptb+2) << "," << *(hptb+3) << ")\n";
        //cout << "g["<<j<<"] = (" << *(gptf-2) << "," << *(gptf-1) << "), g["<<N-j<<"] = (" << *(gptb+2) << "," << *(gptb+3) << ")\n";
    }
    //cout << "N2 = " << N2 << endl;

    *hptf = cos(PI*(N2*N2*omega));
    *(hptf+1) = -sin(PI*(N2*N2*omega));

    //cout << "h["<<N2<<"] = (" << *hptf << "," << *(hptf+1) << ")\n";
    //cout << "g["<<N2<<"] = (" << *gptf << "," << *(gptf+1) << ")\n";

    //for (std::size_t j=0; j!=N; ++j) 
    //    cout << "h["<<j<<"] = ("<< h[2*j] <<","<< h[2*j+1] <<")\n";
    //for (std::size_t j=0; j!=N; ++j) 
    //    cout << "g["<<j<<"] = ("<< g[2*j] <<","<< g[2*j+1] <<")\n";
    

    //Do convolution using Fourier transform

    _FOURIER_TRANSFORM_POWER_OF_TWO_((ValueType*)g.data(), N, 1);
    _FOURIER_TRANSFORM_POWER_OF_TWO_((ScalarType*)h.data(), N, 1);

    for (std::size_t j=0; j!=N; ++j) {
        ValueType mtemp = g[2*j];
        g[2*j] = h[2*j]*mtemp - h[2*j+1]*g[2*j+1];
        g[2*j+1] = h[2*j]*g[2*j+1] + h[2*j+1]*mtemp;
    }

    _FOURIER_TRANSFORM_POWER_OF_TWO_((ValueType*)g.data(), N, -1);

    for(auto& vi: g)
        vi /= (ScalarType) N; 

    //for (std::size_t j=0; j!=N; ++j) 
    //    cout << "h["<<j<<"] = ("<< h[2*j] <<","<< h[2*j+1] <<")\n";

    //--------------------------------------------------------------------------------------//

    // ----- Multiply by weight factors ------
    ScalarType rtemp, itemp;// , mtemp;
    optf = output;
    gptf = g.begin();
    gptb = g.end()-2;

    // First element has weight 1
    *optf = *gptf; // Real part
    ++optf; ++gptf;
    *optf = *gptf; // Imaginary part
    ++optf; ++gptf;
    
    for (std::size_t k=1; k!=m2; ++k) {
        rtemp = first_element(*optf);     //cos(PI*(k*k*omega));
        itemp = first_element(*(optf+1)); //sin(PI*(k*k*omega));

        *optf = rtemp*(*gptf) - itemp*(*(gptf+1));
        *(optf+1) = rtemp*(*(gptf+1)) + itemp*(*gptf);
        optf += 2; gptf += 2;

        *optb = rtemp*(*gptb) - itemp*(*(gptb+1));
        *(optb+1) = rtemp*(*(gptb+1)) + itemp*(*gptb);
        optb -= 2; gptb -= 2;
    }

    if (m_is_even) {
        rtemp = cos(PI*(m2*m2*omega));
        itemp = sin(PI*(m2*m2*omega));
        *optb = rtemp*(*gptb) - itemp*(*(gptb+1));
        *(optb+1) = rtemp*(*(gptb+1)) + itemp*(*gptb);
    }

    //print_ptr = output;
    //for (std::size_t j=0; j!=n-1; ++j, print_ptr+=2)
    //    std::cout << "output["<<j<<"] = ("<<*print_ptr<<","<<*(print_ptr+1)<<")\n";

    //std::cout << "_FQFT_POWER_OF_TWO_PLUS_ONE_ finished" << std::endl;
}



#endif
