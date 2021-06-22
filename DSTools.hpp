#ifndef DSTOOLS
#define DSTOOLS

#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <vector>
#include <string>
//#include <cassert>
#include<type_traits>
#include "WBtemplate.hpp"
#include "LinAlgTools.hpp"


//--------------------------------------------------------------
//-------- Print functions for orbits and angle functions ------
//--------------------------------------------------------------

template<typename VecArray>
void printOrbit(const VecArray& Orbit) {
    std::ofstream att;
    att.open("orbit.txt");
    for (auto u : Orbit) {
        auto ui = u.cbegin();
        auto ue = u.cend()-1;
        for (; ui != ue; ++ui)
            att << *ui << " ";
        att << *ui << "\n";
    }
}

template<typename ScalArray, typename VecArray>
void printDelta(const ScalArray& Delta, const VecArray& Orbit) {
    std::ofstream lift;
    lift.open("lift.txt");
    auto di = Delta.cbegin();
    auto de = Delta.cend();
    auto u = Orbit.cbegin();
    for (;di!=de;++di,++u) {
        for (auto const ui : *u)
            lift << ui << " ";
        lift << *di << "\n";
    }
}


// ----------- Error handling -----------------------


class Failure {
    std::string* s_ptr;
    public:
    Failure() : s_ptr(nullptr) {}
    Failure(std::string _s) : s_ptr(new std::string(_s)) {}
    ~Failure() {
        if (s_ptr!=nullptr)
            delete s_ptr;
    }
    Failure(Failure&& f) {
        s_ptr = f.s_ptr;
        f.s_ptr = nullptr;
    }
    Failure(const Failure& f) : s_ptr(f.s_ptr ? new std::string(*(f.s_ptr)) : nullptr) {}
    friend void swap(Failure& left, Failure& right) {
        std::swap(left.s_ptr,right.s_ptr);
    }
    Failure& operator=(Failure other) {
        swap(*this,other);
        return *this;
    }

    Failure& append(std::string _s) {
        if (s_ptr == nullptr)
            s_ptr = new std::string(_s);
        else
            s_ptr->append(_s);
        return *this;
    }

    Failure& reset() {
        if (s_ptr != nullptr) {
            delete s_ptr;
            s_ptr = nullptr;
        }
        return *this;
    }

    operator bool() const {
        return bool(s_ptr);
    }

    friend std::ostream& operator<<(std::ostream& stream, Failure& fail) {
        if (fail.s_ptr == nullptr)
            return stream << "No failure";
        else
            return stream << *(fail.s_ptr);
    }
};

template<typename ValueType>
std::string to_string_precision(const ValueType& x) {
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<ValueType>::digits10);
    stream << std::scientific << x;
    return stream.str();
}

namespace test {
    enum class Failure {
        no_failure = 0,
        orbit_escaped,           // In RotNumContour::find_attarctor
        attractor_accuracy,      // In RotNumContour::find_attractor
        angle_continuation,      // In CalcDelta::operator()
        false_position,          // In CaclDelta::false_position
        lost_position            // In CaclDelta::false_position (Should never happen!)
    };

    inline std::ostream& operator<<(std::ostream& stream, Failure fail) {
        switch(fail) {
            case Failure::no_failure:
                return stream << "No failure";
            case Failure::orbit_escaped:
                return stream << "Failure: Orbit escaped set bounds";
            case Failure::attractor_accuracy:
                return stream << "Failure: Unable to find attractor to desired accurcay";
            case Failure::angle_continuation:
                return stream << "Failure: Unable to correct all angles";
            case Failure::false_position:
                return stream << "Failure: Unable to enclose solution to accuracy";
            case Failure::lost_position:
                return stream << "Failure: Solution has fallen outside of enclosure";
        }
    }

}

template<typename ValueType>
struct weights {
    ValueType operator()(ValueType t) {
        using std::exp;
        return exp(ValueType(1)/(t*(t-ValueType(1))));
    }
};



//------------------------  DeltaCircle --------------------------

template<typename ValueType, typename Function>
    class DeltaCircle {
        private:
            Function map;
        public:
            typedef ValueType value_type;

            template<typename ScalArray, typename Parameters>
                Failure operator()(ScalArray& Delta,
                        ScalArray& Orbit,
                        const std::size_t N,
                        const Parameters M) {

                    ValueType u = Orbit[0], u_next;

                    Orbit.resize(N);
                    Delta.resize(N-1);

                    using std::floor;

                    for (std::size_t i=0; i!=N-1; ++i) {
                        u_next = map(u,M);
                        Orbit[i] = u;
                        Delta[i] = u_next - u;
                        u = u_next - floor(u_next);
                    }
                    Orbit[N-1] = u;

                    return Failure();
                }

            template<typename ScalArray, typename Parameters>
                Failure operator()(ScalArray&,
                        ScalArray& OrbitBase,
                        ScalArray& DeltaTarget,
                        ScalArray& OrbitTarget,
                        const Parameters M) {

                    return (*this)(M,DeltaTarget,OrbitTarget,OrbitBase.size()); // Might be wrong?
                }
    };

//---------------DeltaCircleAttractor-------------------------

template<typename ValueType, typename Vector, typename Function>
class DeltaCircleAttractor {
    static_assert(std::is_same<ValueType,
            typename Vector::value_type>::value,
            "For CaclDelta<ValueType,Vector>: Vector::value_type must be ValueType");
    public:
    typedef Vector value_type;

    private:
    // Function describing the family of dynamical systems
    Function map;

    // Parameters
    ValueType orbittol;
    std::size_t ContIter;
    ValueType escbound;
    std::size_t OrbitIter;

    template<typename VecArray, typename FunctionParameters>
        Failure find_attractor(VecArray& Orbit,
                const std::size_t N,
                const FunctionParameters MP);

    public:

    DeltaCircleAttractor(ValueType _orbittol, std::size_t _ContIter,
            ValueType _escbound, std::size_t _OrbitIter)
        : orbittol(_orbittol), ContIter(_ContIter),
        escbound(_escbound), OrbitIter(_OrbitIter) { }

    template<typename ScalArray, typename VecArray,
        typename FunctionParameters>
            Failure operator()(ScalArray& Delta,
                    VecArray& Orbit,
                    const FunctionParameters M) {
                Failure fail;

                //std::cout << "Orbit[0] = " << Orbit[0] << std::endl;

                fail = find_attractor(Orbit,Orbit.size(),M);
                if ( fail )
                    return fail.append(": In DeltaCircleAttractor()");
                Delta.resize(Orbit.size()-1);
                for(std::size_t k=0; k != Delta.size(); ++k) {
                    Delta[k] = Orbit[k+1][0]-Orbit[k][0];
                    if(Delta[k] < ValueType(0)) Delta[k] += ValueType(1);
                }

                return fail;
            }
};

template<typename ValueType, typename Vector, typename Function>
template<typename VecArray, typename FunctionParameters>
Failure DeltaCircleAttractor<ValueType,Vector,Function>::find_attractor(VecArray& Orbit,
        const std::size_t N,
        const FunctionParameters MP) {

    static_assert(std::is_same< Vector, typename VecArray::value_type >::value,
            "VecArray::value_type must be Vector");

    Failure fail;

    Vector u = Orbit[0], u0;
    ValueType err = ValueType(2)*orbittol;

    //std::cout << "u = " << u << std::endl; //----------------------
    //std::cout << "MP = (" << MP.mu << "," << MP.epsilon << "," << MP.lambda << ")" << std::endl;

    for (std::size_t j=0; j!=ContIter && norm_1(u) < escbound; ++j)
        u = map(u,MP);

    //std::cout << "u = " << u << std::endl; //----------------------
    //std::cout << "norm_1(u) = " << norm(u) << " and escbound = " << escbound << std::endl; //-------------

    for (std::size_t i=0; i!=OrbitIter && err > orbittol && norm_1(u) < escbound; ++i) {
        u0 = u;
        for (std::size_t j=0; j!=ContIter && err > orbittol && norm_1(u) < escbound; ++j) {
            //std::cout << "u = " << u << "\n";
            u = map(u,MP);
            err = norm_1(u-u0);
        }
    }
    if (!(norm_1(u) < escbound)) {
        fail.append("Failure: Orbit escaped set bounds in find_attractor");
    }
    else if(err > orbittol) {
        fail.append("Failure: Unable to find attractor to desired accuracy");
    }

    Orbit.resize(N);
    Orbit[0] = u;
    for (auto ui=Orbit.begin()+1; ui!=Orbit.end(); ++ui)
        *ui = map(*(ui-1),MP);

    return fail;
}

//------------RotNumCircleAttractor-------------------------------

template<typename ValueType, typename Vector, typename Function>
class RotNumCircleAttractor {
    static_assert(std::is_same<ValueType,
            typename Vector::value_type>::value,
            "For RotNumCircleAttractor<ValueType,Vector,Function>: Vector::value_type must be ValueType");

    DeltaCircleAttractor<ValueType, Vector, Function> CalcDelta;
    WBcoeff< ValueType, weights<ValueType> > coeff;
    std::vector<ValueType> Delta;
    public:

    typedef ValueType value_type;

    class state {
        std::vector<Vector> Orbit;
        public:
        state(std::size_t N) : Orbit(N) {}

        std::vector<Vector>& getOrbit() { return Orbit; }
        const std::vector<Vector>& getOrbit() const { return Orbit; }
        void resize(std::size_t N) { Orbit.resize(N); }
        friend void swap(state r, state l) {
            using std::swap;
            swap(r.Orbit,l.Orbit);
        }
    };

    RotNumCircleAttractor(ValueType orbittol, std::size_t ContIter,
            ValueType escbound, std::size_t OrbitIter) :
        CalcDelta(orbittol, ContIter, escbound, OrbitIter), coeff(), Delta() {}

    template<typename FunctionParameters>
        Failure operator()(ValueType& RotNum, state& State, FunctionParameters FP) {
            Failure fail = CalcDelta(Delta,State.getOrbit(),FP);
            RotNum = WBsum(coeff, Delta.begin(), Delta.size());
            return fail;
        }

    template<typename FunctionParameters>
        Failure operator()(ValueType& RotNum, const state& Base, state& Target, FunctionParameters FP) {
            Target.getOrbit()[0] = Base.getOrbit()[0];
            //std::cout << "Base.getOrbit()[0] = " << Base.getOrbit()[0] << std::endl;
            Failure fail = CalcDelta(Delta,Target.getOrbit(),FP);
            RotNum = WBsum(coeff, Delta.begin(), Delta.size());
            return fail;
        }
};

//---------------- Vectors for projections used in DeltaAttractor -------------------

template<typename Vector>
class ProjVectors{
    private:
        Vector vec_c;
        Vector vec_u;
        Vector vec_v;
    public:
        ProjVectors(Vector _c, Vector _u, Vector _v) :
            vec_c{std::forward<Vector>(_c)},
            vec_u{std::forward<Vector>(_u)},
            vec_v{std::forward<Vector>(_v)}
        {}
        Vector& c() { return vec_c; }
        Vector& u() { return vec_u; }
        Vector& v() { return vec_v; }
        const Vector& c() const { return vec_c; }
        const Vector& u() const { return vec_u; }
        const Vector& v() const { return vec_v; }
};

//------------------------  DeltaAttractor ------------------------------------------------
// Takes a orbit on an attracting circle and calculates angle coordinates for the points on
// the orbit. These angles are used by RotNumAttractor to calculate the rotation number on
// the attractor.

template<typename ValueType, typename Vector, typename Function>
class DeltaAttractor {
    static_assert(std::is_same<ValueType,
            typename Vector::value_type>::value,
            "For CaclDelta<ValueType,Vector>: Vector::value_type must be ValueType");
    public:
    typedef Vector value_type;

    private:
    // Function describing the family of dynamical systems
    Function map;

    // Parameters for numerical method
    std::size_t Lag;
    std::size_t MaxCorr;
    ValueType conttol;
    ProjVectors<Vector> PV;
    ValueType orbittol;
    std::size_t ContIter;
    ValueType escbound;
    std::size_t OrbitIter;

    // Internal working memory
    std::vector<int> checked;
    std::vector<ValueType> u;
    std::vector<ValueType> v;


    template<typename VecArray, typename FunctionParameters>
        Failure find_attractor(VecArray& Orbit,
                const std::size_t N,
                const FunctionParameters MP);

    template<typename ScalArray, typename VecArray>
        void getAngles(ScalArray& Delta,
                const VecArray& Orbit);

    template<typename ScalArray, typename VecArray>
        void continuation(ScalArray& Delta,
                const VecArray& Orbit,
                std::size_t& num_checked);

    template<typename ScalArray, typename VecArray>
        void adaption(const ScalArray& DeltaBase,
                const VecArray& OrbitBase,
                ScalArray& DeltaTarget,
                VecArray& OrbitTarget,
                std::size_t& num_checked);

    public:

    DeltaAttractor(std::size_t _MaxCorr,
            ValueType _conttol, ProjVectors<Vector> _PV,
            ValueType _orbittol, std::size_t _ContIter,
            ValueType _escbound, std::size_t _OrbitIter)
        : Lag(0), MaxCorr(_MaxCorr),
        conttol(_conttol), PV(_PV),
        orbittol(_orbittol), ContIter(_ContIter),
        escbound(_escbound), OrbitIter(_OrbitIter) { }

    template<typename ScalArray, typename VecArray,
        typename FunctionParameters>
            Failure operator()(ScalArray& Delta,
                    VecArray& Orbit,
                    std::size_t N,
                    const std::size_t _Lag,
                    const FunctionParameters M,
                    WBcoeff< ValueType, weights<ValueType> > coeff);

    template<typename ScalArray, typename VecArray,
        typename FunctionParameters>
            Failure operator()(const ScalArray& DeltaBase,
                    const VecArray& OrbitBase,
                    ScalArray& DeltaTarget,
                    VecArray& OrbitTarget,
                    const std::size_t _Lag,
                    const FunctionParameters M);
};

template<typename ValueType, typename Vector, typename Function>
template<typename ScalArray, typename VecArray, typename FunctionParameters>
Failure DeltaAttractor<ValueType,Vector,Function>::operator()(ScalArray& Delta,
        VecArray& Orbit,
        std::size_t N,
        const std::size_t _Lag,
        const FunctionParameters M,
        WBcoeff< ValueType, weights<ValueType> > coeff) {

    static_assert(std::is_same< ValueType, typename ScalArray::value_type >::value,
            "ScalArray::value_type must be ValueType");
    static_assert(std::is_same< Vector, typename VecArray::value_type >::value,
            "VecArray::value_type must be Vector");

    Lag = _Lag;
    Failure fail;

    fail = find_attractor(Orbit,N,M);
    if ( fail )
        return fail.append(": In DeltaAttractor()");

    getAngles(Delta, Orbit);

    // Correct angles in Delta
    checked.resize(Delta.size());
    for (auto&& ch : checked)
        ch = 0;
    checked[0] = 1;
    std::size_t num_checked = 1;
    continuation(Delta, Orbit, num_checked);

    if (num_checked != checked.size()) {
        fail.append("Failure: Unable to correct all angles");
        fail.append(", with num_checked = ").append(std::to_string(num_checked));
        fail.append(" and checked.size() = ").append(std::to_string(checked.size()));
        fail.append(": In DeltaAttractor()");
        //std::cout << "num_checked = " << num_checked << std::endl;
        //std::cout << "checked.size() = " << checked.size() << std::endl;
        // for (size_t i=0; i!=checked.size(); ++i)
        //   if (!checked[i])
        // 	std::cout << "i = " << i << std::endl;
    }

    ValueType rho = WBsum(coeff,Delta.cbegin(),Delta.size());

    using std::floor;
    int r = floor(rho);
    if (r!=0)
        for (auto&& D : Delta)
            D -= r;

    return fail;
}

template<typename ValueType, typename Vector, typename Function>
template<typename ScalArray, typename VecArray, typename FunctionParameters>
Failure DeltaAttractor<ValueType,Vector,Function>::operator()(const ScalArray& DeltaBase,
        const VecArray& OrbitBase,
        ScalArray& DeltaTarget,
        VecArray& OrbitTarget,
        const std::size_t _Lag,
        const FunctionParameters M) {

    static_assert(std::is_same< ValueType, typename ScalArray::value_type >::value,
            "ScalArray::value_type must be ValueType");
    static_assert(std::is_same< Vector, typename VecArray::value_type >::value,
            "VecArray::value_type must be Vector");

    Lag = _Lag;
    Failure fail;

    OrbitTarget.front() = OrbitBase.back();
    fail = find_attractor(OrbitTarget,OrbitBase.size(),M);
    if ( fail )
        return fail.append(": In DeltaAttractor() (adapted)");

    getAngles(DeltaTarget, OrbitTarget);

    // Correct angles in DeltaTarget
    checked.resize(DeltaTarget.size());
    for (auto&& ch : checked)
        ch = 0;
    std::size_t num_checked = 0;

    adaption(DeltaBase, OrbitBase, DeltaTarget, OrbitTarget, num_checked);
    if (num_checked == 0) {
        checked[0] = 1;
        num_checked = 1;
    }

    continuation(DeltaTarget, OrbitTarget, num_checked);

    if (num_checked != checked.size()) {
        fail.append("Failure: Unable to correct all angles");
        fail.append(", with num_checked = ").append(std::to_string(num_checked));
        fail.append(" and checked.size() = ").append(std::to_string(checked.size()));
        fail.append(": In DeltaAttractor() (adapted)");
    }

    return fail;
}

template<typename ValueType, typename Vector, typename Function>
template<typename VecArray, typename FunctionParameters>
Failure DeltaAttractor<ValueType,Vector,Function>::find_attractor(VecArray& Orbit,
        const std::size_t N,
        const FunctionParameters MP) {

    static_assert(std::is_same< Vector, typename VecArray::value_type >::value,
            "VecArray::value_type must be Vector");

    Failure fail;

    Vector u = Orbit[0], u0;
    ValueType err = ValueType(2)*orbittol;

    for (std::size_t j=0; j!=ContIter && norm(u) < escbound; ++j)
        u = map(u,MP);

    for (std::size_t i=0; i!=OrbitIter && err > orbittol && norm(u) < escbound; ++i) {
        u0 = u;
        for (std::size_t j=0; j!=ContIter && err > orbittol && norm(u) < escbound; ++j) {
            u = map(u,MP);
            err = norm(u-u0);
        }
    }
    if (!(norm(u) < escbound)) {
        fail.append("Failure: Orbit escaped set bounds");
    }
    else if(err > orbittol) {
        fail.append("Failure: Unable to find attractor to desired accuracy");
    }

    Orbit.resize(N);
    Orbit[0] = u;
    for (auto ui=Orbit.begin()+1; ui!=Orbit.end(); ++ui)
        *ui = map(*(ui-1),MP);

    return fail;
}

template<typename ValueType, typename Vector, typename Function>
template<typename ScalArray, typename VecArray>
void DeltaAttractor<ValueType,Vector,Function>::getAngles(ScalArray& Delta,
        const VecArray& Orbit) {
    static_assert(std::is_same<ValueType,
            typename ScalArray::value_type>::value,
            "ScalArray::value_type must be ValueType");
    static_assert(std::is_same<Vector,
            typename VecArray::value_type>::value,
            "VecArray::value_type must be Vector");

    std::size_t ND = Orbit.size() - std::max(Lag,std::size_t(1));

    u.resize(ND+1), v.resize(ND+1);

    typename std::vector<ValueType>::iterator ui = u.begin(), vi = v.begin();
    typename VecArray::const_iterator oi = Orbit.begin();

    for(; ui!=u.end(); ++oi, ++ui, ++vi) {
        *ui = prod(*oi-PV.c(),PV.u());
        *vi = prod(*oi-PV.c(),PV.v());
    }

    Delta.resize(ND);

    ui = u.begin(); vi = v.begin();
    const ValueType PI2 = M_PI*2;
    using std::atan2;
    for (auto&& D : Delta)
        D = (atan2(*ui,*vi) - atan2(*(++ui),*(++vi)))/PI2;
}

template<typename ValueType, typename Vector, typename Function>
template<typename ScalArray, typename VecArray>
void DeltaAttractor<ValueType,Vector,Function>::continuation(ScalArray& Delta,
        const VecArray& Orbit,
        std::size_t& num_checked) {
    static_assert(std::is_same<ValueType,
            typename ScalArray::value_type>::value,
            "ScalArray::value_type must be ValueType");
    static_assert(std::is_same<Vector,
            typename VecArray::value_type>::value,
            "VecArray::value_type must be Vector");

    // std::cout << "checked.size() = " << checked.size();
    // std::cout << ", Delta.size() = " << Delta.size();
    // std::cout << ", Orbit.size() = " << Orbit.size() << std::endl;

    using std::round;

    ValueType dist = 0.0;
    std::size_t iter = 0, N = checked.size();
    typename ScalArray::iterator Delta_i;
    typename ScalArray::iterator Delta_j;
    typename VecArray::const_iterator Orbit_i;
    typename VecArray::const_iterator Orbit_j;
    while (num_checked != N && iter < MaxCorr) {

        // std::cout << " --- num_checked = " << num_checked << std::endl;

        Delta_i = Delta.begin();
        Orbit_i = Orbit.begin();
        for (auto checked_i=checked.begin(); checked_i!=checked.end();
                ++checked_i, ++Delta_i, ++Orbit_i) {
            if (*checked_i == 1) {
                Delta_j = Delta.begin();
                Orbit_j = Orbit.begin();
                for (auto checked_j=checked.begin(); checked_j!=checked.end();
                        ++checked_j, ++Delta_j, ++Orbit_j) {
                    if (*checked_j == 0) {
                        dist = norm((*Orbit_i)-(*Orbit_j));
                        for (std::size_t l=1; l!=Lag+1; ++l)
                            dist +=  norm(*(Orbit_i+l)-*(Orbit_j+l));
                        dist /= Lag+1;
                        if (dist < conttol) {
                            *Delta_j -= round((*Delta_j)-(*Delta_i));
                            *checked_j = 1; ++num_checked;
                            // std::cout << "num_checked = " << num_checked;
                            // std::cout << ", i = " << checked_i-checked.begin();
                            // std::cout << ", j= " << checked_j-checked.begin() << std::endl;
                        }
                    }
                }
            }
        }
        ++iter;
    }
}

template<typename ValueType, typename Vector, typename Function>
template<typename ScalArray, typename VecArray>
void DeltaAttractor<ValueType,Vector,Function>::adaption(const ScalArray& DeltaBase,
        const VecArray& OrbitBase,
        ScalArray& DeltaTarget,
        VecArray& OrbitTarget,
        std::size_t& num_checked) {
    static_assert(std::is_same<ValueType,
            typename ScalArray::value_type>::value,
            "ScalArray::value_type must be ValueType");
    static_assert(std::is_same<Vector,
            typename VecArray::value_type>::value,
            "VecArray::value_type must be Vector");

    ValueType dist;
    std::size_t N = checked.size();
    for (std::size_t j=0; j!=N; ++j)
        if (!checked[j])
            for (std::size_t i=0; i!=N; ++i) {
                dist = norm(OrbitBase[i]-OrbitTarget[j]);
                for (std::size_t l=1; l!=Lag+1; ++l)
                    dist += norm(OrbitBase[i+l]-OrbitTarget[j+l]);
                dist /= Lag+1;
                if (dist < conttol) {
                    DeltaTarget[j] -= std::round(DeltaTarget[j]-DeltaBase[i]);
                    checked[j] = 1; ++num_checked;
                    break;
                }
            }
}


//------------------------  RotNumAttractor ------------------------
// Calculates rotation numer for circle type attractors. Stores orbit data in States
// that can be saved and reused between calculation for efficiency.

template<typename ValueType, typename Vector, typename Function>
class RotNumAttractor {
    static_assert(std::is_same<ValueType,
            typename Vector::value_type>::value,
            "For RotNumAttractor<ValueType,Vector,CalcDelta>: Vector::value_type must be ValueType");

    DeltaAttractor<ValueType, Vector, Function> CalcDelta;
    WBcoeff< ValueType, weights<ValueType> > coeff;
    public:

    typedef ValueType value_type;

    class state {
        std::vector<ValueType> Delta;
        std::vector<Vector> Orbit;
        std::size_t Lag;
        public:
        state(std::size_t N, std::size_t _Lag) : Delta(N - std::max(_Lag,std::size_t(1))), Orbit(N), Lag(_Lag) {}

        std::vector<ValueType>& getDelta() { return Delta; }
        std::vector<Vector>& getOrbit() { return Orbit; }
        const std::vector<ValueType>& getDelta() const { return Delta; }
        const std::vector<Vector>& getOrbit() const { return Orbit; }

        std::size_t getLag() const { return Lag;}
        void setLag(std::size_t L) {
            Delta.resize(Orbit.size() - std::max(L,std::size_t(1)));
            Lag = L;
        }
        void resize(std::size_t N) {
            Delta.resize(N - std::max(Lag,std::size_t(1)));
            Orbit.resize(N);
        }
    };

    RotNumAttractor(std::size_t MaxCorr,
            ValueType conttol, ProjVectors<Vector> PV,
            ValueType orbittol, std::size_t ContIter,
            ValueType escbound, std::size_t OrbitIter) :
        CalcDelta(MaxCorr,conttol,PV, orbittol,ContIter,escbound,OrbitIter), coeff() {}

    // state makestate(std::size_t N, std::size_t Lag) {
    //   return std::move(state(N,Lag));
    // }

    template<typename FunctionParameters>
        Failure operator()(ValueType& RotNum, state& State, FunctionParameters FP) {
            Failure fail = CalcDelta(State.getDelta(),State.getOrbit(),State.getOrbit().size(),State.getLag(),FP,coeff);
            RotNum = WBsum(coeff, State.getDelta().begin(), State.getDelta().size());
            return fail;
        }

    template<typename FunctionParameters>
        Failure operator()(ValueType& RotNum, const state& Base, state& Target, FunctionParameters FP) {
            Failure fail = CalcDelta(Base.getDelta(),Base.getOrbit(),Target.getDelta(),Target.getOrbit(),Base.getLag(),FP);
            RotNum = WBsum(coeff, Target.getDelta().begin(), Target.getDelta().size());
            return fail;
        }

    WBcoeff< ValueType, weights<ValueType> >& getcoeff() {
        return coeff;
    }
};

template<typename ValueType, typename Vector, typename Function>
void swap(typename RotNumAttractor<ValueType,Vector,Function>::state& left, typename RotNumAttractor<ValueType,Vector,Function>::state& right) {
    std::swap(left.getDelta(),right.getDelta());
    std::swap(left.getOrbit(),right.getOrbit());
    std::size_t Lag = left.getLag();
    left.setLag(right.getLag());
    right.setLag(Lag);
}


//------------------------  FunctionCountor --------------------------
// Basic tools for finding points for which a given function has certain values.
// Searches for values along lines or circles, and can be used to trace countors.
// Assumes the function has internal states that are saved and reused between calls.

template<typename ValueType>
class FunctionContour {
    public:
        class iterator;

    private:
        // Parameters for numerical method
        ValueType stepsize;
        ValueType paramtol;
        std::size_t MaxBis;
        ValueType rhotol;
        std::size_t NumSearch;
        ValueType searchrad;

    public:
        FunctionContour(ValueType _stepsize,
                ValueType _paramtol,
                std::size_t _MaxBis,
                ValueType _rhotol,
                std::size_t _NumSearch,
                ValueType _searchrad) :
            stepsize(_stepsize), paramtol(_paramtol),
            MaxBis(_MaxBis), rhotol(_rhotol),
            NumSearch(_NumSearch), searchrad(_searchrad) {}

        void setsearchrad(ValueType _searchrad) {
            searchrad = _searchrad;
        }

        template<typename Function, typename Parameters>
            Failure search_line(Parameters& ans,
                    const ValueType rho,
                    Function& function,
                    typename Function::state& Front,
                    typename Function::state& Back,
                    Parameters begin,
                    Parameters end);

        template<typename Function, typename Parameters>
            Failure search_circle(Parameters& ans,
                    const ValueType rho,
                    Function& function,
                    const typename Function::state& Base,
                    typename Function::state& Target1,
                    typename Function::state& Target2,
                    const Parameters Pbase,
                    const Parameters Pprev);

    private:

        template<typename Function, typename FunctionParameters>
            Failure false_position(FunctionParameters& Mmid,
                    const ValueType rho,
                    Function& function,
                    ValueType rholeft,
                    FunctionParameters Mleft,
                    ValueType rhoright,
                    FunctionParameters Mright,
                    typename Function::state& Target,
                    typename Function::state& Base);

        template<typename FunctionParameters>
            FunctionParameters mid(FunctionParameters left, FunctionParameters& right) {
                left[0] = ( left[0] + right[0] )/ValueType(2);
                left[1] = ( left[1] + right[1] )/ValueType(2);
                return left;
            }

        template<typename FunctionParameters>
            FunctionParameters false_mid(FunctionParameters left, FunctionParameters& right,
                    ValueType Fleft, ValueType Fright) {
                left[0] = ( Fright*left[0] - Fleft*right[0] )/(Fright-Fleft);
                left[1] = ( Fright*left[1] - Fleft*right[1] )/(Fright-Fleft);
                return left;
            }

        template<typename FunctionParameters>
            ValueType forward(FunctionParameters Mprev,
                    FunctionParameters Mbase,
                    FunctionParameters Mcand) {
                Mprev[0] -= Mbase[0]; Mprev[1] -= Mbase[1];
                Mcand[0] -= Mbase[0]; Mcand[1] -= Mbase[1];
                ValueType ans = Mprev[0]*Mcand[0] + Mprev[1]*Mcand[1];
                using std::sqrt;
                ans /= sqrt(Mprev[0]*Mprev[0] + Mprev[1]*Mprev[1]);
                ans /= sqrt(Mcand[0]*Mcand[0] + Mcand[1]*Mcand[1]);
                return -ans;
            }
};


template<typename ValueType>
template<typename Function, typename Parameters>
Failure FunctionContour<ValueType>::search_line(Parameters& ans,
        const ValueType rho,
        Function& function,
        typename Function::state& Front,
        typename Function::state& Back,
        Parameters begin,
        Parameters end) {
    static_assert(std::is_same<ValueType,
            typename Function::value_type>::value,
            "Function::value_type must be ValueType");

    using std::ceil;
    using std::sqrt;
    using std::swap;

    Failure fail;

    ValueType rhofront(0), rhoback(0);
    Parameters Pfront = begin, Pback;

    fail = function(rhofront, Front, Pfront);

    if ( fail )
        return fail.append(": In FunctionContour::search_line");

    ValueType lenP1 = end[0]-begin[0], lenP2 = end[1]-begin[1];
    ValueType lenP = sqrt(lenP1*lenP1+lenP2*lenP2)/stepsize;
    std::size_t ParNum = ceil(lenP);

    // std::cout << "lenP1 = " << lenP1 << ", lenP2 = " << lenP2 << std::endl;
    // std::cout << "lenP = " << lenP << ", ParNum = " << ParNum << std::endl;
    // std::cout << "rho = " << rho << std::endl;


    for (std::size_t m=0; m!=ParNum; ++m) {
        Pback = Pfront;
        Pfront[0] = begin[0]+lenP1*m/lenP;
        Pfront[1] = begin[1]+lenP2*m/lenP;

        swap(rhofront,rhoback);
        swap(Front,Back);

        fail = function(rhofront,Back,Front,Pfront);
        if ( fail ) {
            fail.append(": In FunctionContour::search_line with m = ").append(std::to_string(m));
            fail.append(" and ParNum = ").append(std::to_string(ParNum));
            return fail;
            break;
        }

        if ((rhoback <= rho && rho <= rhofront)||(rhofront <= rho && rho <= rhoback))
            break;

        if ( m+1 == ParNum )
            fail.append("Failure: Unable to find value on line");
    }

    if ( fail )
        return fail;

    return false_position(ans,rho,function,rhoback,Pback,rhofront,Pfront,Back,Front);
}

template<typename ValueType>
template<typename Function, typename Parameters>
Failure FunctionContour<ValueType>::false_position(Parameters& Mmid,
        const ValueType rho,
        Function& function,
        ValueType rholeft,
        Parameters Mleft,
        ValueType rhoright,
        Parameters Mright,
        typename Function::state& Target,
        typename Function::state& Base) {
    static_assert(std::is_same<ValueType,
            typename Function::value_type>::value,
            "Function::value_type must be ValueType");

    using std::sqrt;
    using std::swap;
    using std::abs;

    Failure fail;

    ValueType rhomid = rholeft;
    std::size_t bisiter = 0;
    ValueType Mdist = sqrt(  (Mleft[0]-Mright[0])
            *(Mleft[0]-Mright[0])
            + (Mleft[1]-Mright[1])
            *(Mleft[1]-Mright[1]) );


    while ( Mdist > paramtol && abs(rho-rhomid) > rhotol && bisiter < MaxBis && !fail) {
        Mmid = false_mid(Mleft,Mright,rholeft-rho,rhoright-rho);

        //std::cout << "rholeft-rho = " << rholeft - rho << ", rhoright-rho = " << rhoright-rho << std::endl;
        //std::cout << "Mmid[0] = " << Mmid[0] << ", Mmid[1] = " << Mmid[1] << std::endl;

        fail = function(rhomid, Base, Target, Mmid);

        if ( fail ) {
            fail.append(": In FunctionContour::false_position with bisiter = ").append(std::to_string(bisiter));
        }

        if ((rholeft-rhomid)*(rhoright-rhomid) > ValueType(0)) {
            fail.append("Failure: Unable to proceed in false_position");
            fail.append(" with rholeft - rhomid = ").append(to_string_precision(rholeft-rhomid));
            fail.append(" and rhoright - rhomid = ").append(to_string_precision(rhoright-rhomid));
            fail.append(", and bisiter = ").append(std::to_string(bisiter));
            break;
        }

        if ((rho-rhomid)*(rho-rholeft) >= ValueType(0)) {
            Mleft = Mmid;
            rholeft = rhomid;
        } else if ((rho-rhomid)*(rho-rhoright) >= ValueType(0)) {
            Mright = Mmid;
            rhoright = rhomid;
        } else {
            fail.append("Failure: Solution has fallen outside of enclosure");
            fail.append(" with rholeft - rho = ").append(to_string_precision(rholeft-rho));
            fail.append(" and rhoright - rho = ").append(to_string_precision(rhoright-rho));
            break;
        }

        swap(Target,Base);

        Mdist = sqrt(  (Mleft[0]-Mright[0])
                *(Mleft[0]-Mright[0])
                + (Mleft[1]-Mright[1])
                *(Mleft[1]-Mright[1]) );
        ++bisiter;
    }

    if ( fail )
        return fail.append(": In FunctionContour::false_position");

    //std::cout << "bisiter = " << bisiter << std::endl;

    if ( abs(rho-rhomid) > rhotol ) {
        std::cout << "Where?" << std::endl;
        fail.append("Failure: Unable to enclose solution to accuracy with error = ");
        fail.append(to_string_precision(abs(rho-rhomid)));
        fail.append(", and bisiter = ").append(std::to_string(bisiter));
        fail.append(": In false_position");
    }

    return fail;
}

template<typename ValueType>
template<typename Function, typename Parameters>
Failure FunctionContour<ValueType>::search_circle(Parameters& ans,
        const ValueType rho,
        Function& function,
        const typename Function::state& Base,
        typename Function::state& Target1,
        typename Function::state& Target2,
        const Parameters Pbase,
        const Parameters Pprev) {
    static_assert(std::is_same<ValueType,
            typename Function::value_type>::value,
            "Function::value_type must be ValueType");

    using std::cos;
    using std::sin;
    using std::atan2;
    using std::swap;

    Failure fail;

    ValueType rhoright, rholeft;
    ValueType start_angle = atan2(Pbase[1]-Pprev[1],Pbase[0]-Pprev[0]); //Used to be atan2(Pbase[0]-Pprev[0],Pbase[1]-Pprev[1])
    Parameters Pright = Pbase, Pleft = Pbase;
    Pright[0] = Pbase[0] + searchrad*cos(start_angle - M_PI/NumSearch);
    Pright[1] = Pbase[1] + searchrad*sin(start_angle - M_PI/NumSearch);
    Pleft[0] = Pbase[0] + searchrad*cos(start_angle + M_PI/NumSearch);
    Pleft[1] = Pbase[1] + searchrad*sin(start_angle + M_PI/NumSearch);

    fail = function(rhoright, Base, Target1, Pright);
    if ( fail ) {
        return fail.append(": In search_circle for rhoright ");
    }

    fail = function(rholeft, Base, Target2, Pleft);
    if ( fail ) {
        return fail.append(": In search_circle for rholeft ");
    }

    int direction = 0;
    if ( (rholeft - rho)*(rhoright - rho) > 0 ) {
        if ( (rhoright - rholeft)*(rho - rholeft) > 0) {
            direction = -1;
            swap(Pright,Pleft);
            swap(Target1,Target2);
        } else if  ( (rhoright - rholeft)*(rho - rholeft) < 0) {
            direction = 1;
        } else {
            return fail.append("Failure: Unclear which direction to search: In search_circle");
        }

        std::cout << "start_angle/(2*M_PI) = atan2("<< Pbase[1]-Pprev[1] <<","<<Pbase[0]-Pprev[0]<<")/(2*M_PI) = " << start_angle/(2*M_PI) << "\n"; //--------------------

        for (std::size_t i=0; i!=NumSearch/4+1; ++i) {
            swap(Pright,Pleft);
            swap(Target1,Target2);
            Pleft[0] = Pbase[0] + searchrad*cos(start_angle + direction*(2*M_PI*i)/NumSearch);
            Pleft[1] = Pbase[1] + searchrad*sin(start_angle + direction*(2*M_PI*i)/NumSearch);
            rhoright = rholeft;
            fail = function(rholeft, Base, Target2, Pleft);
            //std::cout << "Pleft = (" << Pleft[0] << ", " << Pleft[1] << ", " << Pleft[2] << ")\n"; // -------------
            //std::cout << "rholeft = " << rholeft << "\n"; // ---------------------
            if ( fail ) {
                return fail.append(": In search_circle with i = ").append(std::to_string(i));
            }
            if ( (rholeft - rho)*(rhoright - rho) < 0 )
                break;
            if ( i == NumSearch/4 )
                return fail.append("Failure: Unable to find value on circle: In search_circle");
        }
    }

    return false_position(ans,rho,function,rholeft,Pleft,rhoright,Pright,Target2,Target1);
}

// template<typename ValueType>
// template<typename Function, typename Parameters>
// Failure FunctionContour<ValueType>::search_circle(Parameters& ans,
// 						  const ValueType rho,
// 						  Function& function,
// 						  typename Function::state& Base,
// 						  typename Function::state& Target,
// 						  Parameters Pbase,
// 						  Parameters Pprev) {
//   static_assert(std::is_same<ValueType,
// 		typename Function::value_type>::value,
// 		"Function::value_type must be ValueType");

//   using std::cos;
//   using std::sin;

//   Failure fail;

//   std::vector<Parameters> Pcircle(NumSearch);
//   std::vector<ValueType> rhocircle(NumSearch);
//   for (std::size_t i=0; i!=NumSearch; ++i) {
//     Pcircle[i] = Pbase;
//     Pcircle[i][0] = Pbase[0] + searchrad*cos((2*M_PI*i)/NumSearch);
//     Pcircle[i][1] = Pbase[1] + searchrad*sin((2*M_PI*i)/NumSearch);

//     fail = function(rhocircle[i], Base, Target, Pcircle[i]);
//     if ( fail ) {
//       return fail.append(": In search_circle with i = ").append(std::to_string(i));
//     }
//   }

//   ValueType Pforward, maxforward = -2.0;
//   std::size_t IndexMax = NumSearch;
//   for (std::size_t i=1; i!=NumSearch; ++i) {
//     if ((rho-rhocircle[i])*(rho-rhocircle[i-1]) < 0) {
//       Pforward = forward(Pprev,Pbase,mid(Pcircle[i],Pcircle[i-1]));
//       //std::cout << "Pforward = " << Pforward << std::endl;
//       if (Pforward > maxforward) {
// 	maxforward = Pforward;
// 	IndexMax = i;
//       }
//     }
//   }

//   if ( IndexMax == NumSearch ) {
//     return fail.append("Failure: Unable to find appropriate value on circle");
//   }

//   return false_position(ans,rho,function,rhocircle[IndexMax-1],Pcircle[IndexMax-1],rhocircle[IndexMax],Pcircle[IndexMax],Target,Base);
// }

template<typename ValueType>
class FunctionContour<ValueType>::iterator {


};


#endif
