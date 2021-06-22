
#include <iostream>
#include <iomanip>
#include "DSTools.hpp"
#include "Interval.hpp"

using Real = my_mpfr;
//using Real = double;
const std::size_t prec = 210;

const std::size_t dim = 2;
using Vec2 = Vec<Real,dim>;

template<typename T>
void init() = delete;

template<>
void init<double>() {
    std::cout << "Using double" << std::endl;
    std::cout << std::setprecision(std::numeric_limits<double>::digits10);
}

template<>
void init<my_mpfr>() {
    //std::size_t prec = 100;
    mpfr_set_default_prec(prec);  
    printf("Using my_pfr with precision %lu (MPFR_PREC_MIN = %d, MPFR_PREC_MAX = %ld)\n",prec,MPFR_PREC_MIN,MPFR_PREC_MAX);
}


using std::floor;
using std::sin;

template<typename ValueType>
class MapParameters{
    public:
        ValueType mu, epsilon, lambda;
        ValueType& operator[](std::size_t i){
            return *((ValueType*)this + i);
        }
        const ValueType& operator[](std::size_t i) const {
            return *((ValueType*)this + i);
        }
};

template<typename ValueType>
inline ValueType phi(ValueType x, const ValueType& y, const MapParameters<ValueType>& param) {
    const ValueType twopi = ValueType(2)*pi<ValueType>();
    const ValueType twopirec = ValueType(1)/twopi;
    x = param.lambda*y + param.mu + param.epsilon*twopirec*sin(twopi*x);
    return x;
}

template<typename ValueType, std::size_t dim>
inline Vec<ValueType, dim> F_p(Vec<ValueType, dim> u, const MapParameters<ValueType>& param) {
    u[0] = phi(u[0],u[1],param);
    u[1] = u[0];
    return u;
}
template<typename ValueType, std::size_t dim>
inline Vec<ValueType, dim> F(Vec<ValueType, dim> u, const MapParameters<ValueType>& param) {
    u[1] = phi(u[0],u[1],param);
    u[0] += u[1];
    u[0] -= floor(u[0]);
    return u;
}
template<typename ValueType>
struct Function {
    Vec<ValueType,dim> operator() (Vec<ValueType,dim>& u, const MapParameters<ValueType>& param) {
        return F(u, param);
    }
};

int main(void) {

    init<Real>();
    
    std::size_t OrbitSize = 1<<17;//8000;
    std::size_t OrbitIter = 10; //
    Real orbittol = 0x1p-5; //
    Real escbound = 30.0; //
    std::size_t ContIter = 10000; //
    Real stepsize = 0x1p-8; //
    Real searchrad = 0x1p-6;//0x1p-10;//0x1p-16;
    std::size_t NumSearch = 128; //32;
    Real paramtol = 0x1p-40; //
    std::size_t MaxBis = 50; //
    Real rhotol = pow(Real(10),-55); //0x1p-55; //-37

    RotNumCircleAttractor<Real,Vec2,Function<Real>> rotnum(orbittol, ContIter, escbound, OrbitIter);
    RotNumCircleAttractor<Real,Vec2,Function<Real>>::state State1(OrbitSize);

    using std::sqrt;
    Real golden_mean = (sqrt(Real(5.0))-Real(1.0))/Real(2.0);  

    //Real lambda = std::string("0.4"), epsilon = std::string("0.98142138274"); //0.979142644237101;
    Real lambda = std::string("0.4");//, epsilon = std::string("0.5");
    MapParameters<Real> MP = {(1.0-lambda)*golden_mean, Real(0.0), lambda};

    std::cout << "Searching for rotation number " << golden_mean << " within error " << rhotol << "\n";
    std::cout << "lambda  = " << lambda << std::endl;

    std::size_t EpsNum = 9;
    std::vector<Real> epsilon_list(EpsNum);
    for (std::size_t i=0; i!=EpsNum; ++i)
        epsilon_list[i] = Real(i+1)/(EpsNum+1);

    Failure fail;

    //std::cout << "Searching for rotation number " << golden_mean << " within error " << rhotol << ",\n for epsilon = " << epsilon << " and lambda = " << lambda << std::endl;

    FunctionContour<Real> funcont(stepsize, paramtol, MaxBis, rhotol, NumSearch, searchrad);

    RotNumCircleAttractor<Real,Vec2,Function<Real>>::state State2(OrbitSize), State3(OrbitSize);

    std::ofstream parameters("parameter_list.txt");

    for (std::size_t i=0; i!=EpsNum; ++i) {
        MapParameters<Real> begin = {0.0, epsilon_list[i], lambda}, end = {1.0, epsilon_list[i], lambda};
        State1.getOrbit()[0] = {0.0, golden_mean};

        std::cout << "i = " << i << ", epsilon = " << epsilon_list[i] << "\n";
        fail = funcont.search_line(MP, golden_mean, rotnum, State1, State2, begin, end);

        std::cout << fail << std::endl;
        if (!fail) std::cout << "mu = " << MP.mu << std::endl;

        parameters << MP.mu << " " << MP.epsilon << " " << MP.lambda << "\n";
    }


    return 0;
}


