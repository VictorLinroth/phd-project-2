#ifndef my_mpfr_H
#define my_mpfr_H
/*
	Version 1.2 multiprecision.
	
	Last modification: 10:30h, 12-06-2020
*/

#include <ostream>
#include <cmath>
#include <gmp.h>
#include <mpfr.h>

//using namespace std; Jesus Christ...

class my_mpfr;

template<typename T>
T pi() = delete;

template<>
float pi<float>() {return (float)M_PI;}

template<>
double pi<double>() {return M_PI;}

template<>
my_mpfr pi<my_mpfr>();

class my_mpfr
{
public:
    mpfr_t x;
    
    // Regular constructors and destructor
    my_mpfr();
    my_mpfr(const my_mpfr&);
    my_mpfr(const unsigned long int);
    my_mpfr(const signed long int);
    my_mpfr(const double);
    my_mpfr(const long double);
    my_mpfr(const int);
    my_mpfr(const char[], const int);
    my_mpfr(const std::string, const int);
    ~my_mpfr();
    
    // Regular assignments
    my_mpfr& operator = (const my_mpfr&);
    my_mpfr& operator = (const unsigned long int);
    my_mpfr& operator = (const signed long int);
    my_mpfr& operator = (const double);
    my_mpfr& operator = (const long double);
    my_mpfr& operator = (const int);
    my_mpfr& operator = (const std::string);

    // Move constructor and move assignment
    my_mpfr(my_mpfr&&);
    my_mpfr& operator=(my_mpfr&&);

    // Conversion operators
    operator unsigned long int() { return mpfr_get_ui(x,GMP_RNDN); }
    operator signed long int() { return mpfr_get_si(x,GMP_RNDN); }

    // Compound assigment operators
    inline my_mpfr operator+= (const my_mpfr&);
    inline my_mpfr operator+= (const unsigned long int);
    inline my_mpfr operator+= (const signed long int);
    inline my_mpfr operator+= (const double);
    inline my_mpfr operator+= (const int);

    inline my_mpfr operator-= (const my_mpfr&);
    inline my_mpfr operator-= (const unsigned long int);
    inline my_mpfr operator-= (const signed long int);
    inline my_mpfr operator-= (const double);
    inline my_mpfr operator-= (const int);

    inline my_mpfr operator*= (const my_mpfr&);
    inline my_mpfr operator*= (const unsigned long int);
    inline my_mpfr operator*= (const signed long int);
    inline my_mpfr operator*= (const double);
    inline my_mpfr operator*= (const int);

    inline my_mpfr operator/= (const my_mpfr&);
    inline my_mpfr operator/= (const unsigned long int);
    inline my_mpfr operator/= (const signed long int);
    inline my_mpfr operator/= (const double);
    inline my_mpfr operator/= (const int);

    // Efficient swap function
    friend void swap(my_mpfr& P, my_mpfr& Q)
    {
        //std::cout << "Calling my_mpfr::swap" << std::endl;
        mpfr_swap(P.x, Q.x);
    }
};

// Arithmetic operations +,-,*,/
inline my_mpfr operator- (const my_mpfr& P);

inline my_mpfr operator+ (const my_mpfr&, const my_mpfr&);
inline my_mpfr operator+ (const unsigned long int, const my_mpfr&);
inline my_mpfr operator+ (const my_mpfr&, const unsigned long int);
inline my_mpfr operator+ (const signed long int, const my_mpfr&);
inline my_mpfr operator+ (const my_mpfr&, const signed long int);
inline my_mpfr operator+ (const double, const my_mpfr&);
inline my_mpfr operator+ (const my_mpfr&, const double);
inline my_mpfr operator+ (const int, const my_mpfr&);
inline my_mpfr operator+ (const my_mpfr&, const int);

inline my_mpfr operator- (const my_mpfr&, const my_mpfr&);
inline my_mpfr operator- (const unsigned long int, const my_mpfr&);
inline my_mpfr operator- (const my_mpfr&, const unsigned long int);
inline my_mpfr operator- (const signed long int, const my_mpfr&);
inline my_mpfr operator- (const my_mpfr&, const signed long int);
inline my_mpfr operator- (const double, const my_mpfr&);
inline my_mpfr operator- (const my_mpfr&, const double);
inline my_mpfr operator- (const int, const my_mpfr&);
inline my_mpfr operator- (const my_mpfr&, const int);

inline my_mpfr operator* (const my_mpfr&, const my_mpfr&);
inline my_mpfr operator* (const unsigned long int, const my_mpfr&);
inline my_mpfr operator* (const my_mpfr&, const unsigned long int);
inline my_mpfr operator* (const signed long int, const my_mpfr&);
inline my_mpfr operator* (const my_mpfr&, const signed long int);
inline my_mpfr operator* (const double, const my_mpfr&);
inline my_mpfr operator* (const my_mpfr&, const double);
inline my_mpfr operator* (const int, const my_mpfr&);
inline my_mpfr operator* (const my_mpfr&, const int);

inline my_mpfr operator/ (const my_mpfr&, const my_mpfr&);
inline my_mpfr operator/ (const unsigned long int, const my_mpfr&);
inline my_mpfr operator/ (const my_mpfr&, const unsigned long int);
inline my_mpfr operator/ (const signed long int, const my_mpfr&);
inline my_mpfr operator/ (const my_mpfr&, const signed long int);
inline my_mpfr operator/ (const double, const my_mpfr&);
inline my_mpfr operator/ (const my_mpfr&, const double);
inline my_mpfr operator/ (const int, const my_mpfr&);
inline my_mpfr operator/ (const my_mpfr&, const int);

// Relational operators >,<,>=,<=
inline int operator== (const my_mpfr& P, const my_mpfr& Q);
inline int operator!= (const my_mpfr& P, const my_mpfr& Q);

inline int operator> (const my_mpfr& P, const my_mpfr& Q);
inline int operator> (const my_mpfr& P, unsigned long int Q);
inline int operator> (unsigned long int P, const my_mpfr& Q);
inline int operator> (const my_mpfr& P, signed long int Q);
inline int operator> (signed long int P, const my_mpfr& Q);
inline int operator> (const my_mpfr& P, double Q);
inline int operator> (double P, const my_mpfr& Q);

inline int operator< (const my_mpfr& P, const my_mpfr& Q);
inline int operator< (const my_mpfr& P, unsigned long int Q);
inline int operator< (unsigned long int P, const my_mpfr& Q);
inline int operator< (const my_mpfr& P, signed long int Q);
inline int operator< (signed long int P, const my_mpfr& Q);
inline int operator< (const my_mpfr& P, double Q);
inline int operator< (double P, const my_mpfr& Q);

inline int operator>= (const my_mpfr& P, const my_mpfr& Q);
inline int operator>= (const my_mpfr& P, unsigned long int Q);
inline int operator>= (unsigned long int P, const my_mpfr& Q);
inline int operator>= (const my_mpfr& P, signed long int Q);
inline int operator>= (signed long int P, const my_mpfr& Q);
inline int operator>= (const my_mpfr& P, double Q);
inline int operator>= (double P, const my_mpfr& Q);

inline int operator<= (const my_mpfr& P, const my_mpfr& Q);
inline int operator<= (const my_mpfr& P, unsigned long int Q);
inline int operator<= (unsigned long int P, const my_mpfr& Q);
inline int operator<= (const my_mpfr& P, signed long int Q);
inline int operator<= (signed long int P, const my_mpfr& Q);
inline int operator<= (const my_mpfr& P, double Q);
inline int operator<= (double P, const my_mpfr& Q);

// Common functions
inline my_mpfr sin(const my_mpfr& P);
inline my_mpfr cos(const my_mpfr& P);
inline my_mpfr sinh(const my_mpfr& P);
inline my_mpfr cosh(const my_mpfr& P);
inline my_mpfr exp(const my_mpfr& P);
inline my_mpfr exp2(const my_mpfr& P);
inline my_mpfr exp10(const my_mpfr& P);
inline my_mpfr log(const my_mpfr& P);
inline my_mpfr log2(const my_mpfr& P);
inline my_mpfr log10(const my_mpfr& P);
inline my_mpfr sqrt(const my_mpfr& P);
inline my_mpfr abs(const my_mpfr& P);
inline my_mpfr pow(const my_mpfr& P, int i);
inline my_mpfr pow(const my_mpfr& P, const my_mpfr& x);
inline my_mpfr max(const my_mpfr& a, const my_mpfr& b);
inline my_mpfr max(const my_mpfr& a, const my_mpfr& b, const my_mpfr& c);
inline my_mpfr max(const my_mpfr& a, const my_mpfr& b, const my_mpfr& c, const my_mpfr& d);
inline my_mpfr round(const my_mpfr& P);
inline my_mpfr ceil(const my_mpfr& P);
inline my_mpfr floor(const my_mpfr& P);
inline int isfinite(const my_mpfr& P);
inline int isinf(const my_mpfr& P);
inline int isnan(const my_mpfr& P);

// Stream operator overload
std::ostream& operator<<(std::ostream& stream, const my_mpfr& op); 

/*-------------------------------------------------------------------------------------------------------------------------------------------*/

// Regular constructors and destructor

my_mpfr::my_mpfr () 
{
    //std::cout << "Calling constructor my_mpfr() for " << this << std::endl;
    mpfr_init(x);
} 

my_mpfr::my_mpfr(const my_mpfr& P)
{
    //std::cout << "Calling copy constructor my_mpfr(my_mpfr&) for " << this << std::endl;
    mpfr_init(x);
    mpfr_set(x, P.x, GMP_RNDN);
}

my_mpfr::my_mpfr (const unsigned long int a) 
{
    mpfr_init(x);
    mpfr_set_ui(x, a, GMP_RNDN);
}

my_mpfr::my_mpfr (const signed long int a) 
{
    mpfr_init(x);
    mpfr_set_si(x, a, GMP_RNDN);
}

my_mpfr::my_mpfr (double a) 
{
    //std::cout << "Calling copy constructor my_mpfr(double("<<a<<")) for " << this << std::endl;
    mpfr_init(x);
    mpfr_set_d(x, a, GMP_RNDN);
}

my_mpfr::my_mpfr (const long double a) 
{
    mpfr_init(x);
    mpfr_set_ld(x, a, GMP_RNDN);
}

my_mpfr::my_mpfr(const int a)
{
    mpfr_init(x);
    mpfr_set_si(x, (signed long int)a, GMP_RNDN);
}

my_mpfr::my_mpfr(const char a[], const int b = 10)
{
    mpfr_init(x);
    mpfr_set_str(x, a, b, GMP_RNDN);
}

my_mpfr::my_mpfr(const std::string a, const int b = 10)
{
    mpfr_init(x);
    mpfr_set_str(x, a.data(), b, GMP_RNDN);
}

my_mpfr::~my_mpfr()
{
    //std::cout << "Calling destructor ~my_mpfr() for " << this << std::endl;
    if (x->_mpfr_d!=0)
        mpfr_clear(x);
}

// Regular assignments

my_mpfr& my_mpfr::operator = (const my_mpfr& P)
{
    //std::cout << "Calling copy assignemnt =(my_mpfr&) for " << this << std::endl;
    if(this != &P)
    {
	mpfr_set(x, P.x, GMP_RNDN);
    }
    return *this;
}

my_mpfr& my_mpfr::operator = (const unsigned long int P)
{
    mpfr_set_ui(x, P, GMP_RNDN);

    return *this;
}

my_mpfr& my_mpfr::operator = (const signed long int P)
{
    mpfr_set_si(x, P, GMP_RNDN);

    return *this;
}

my_mpfr& my_mpfr::operator = (const double P)
{
    mpfr_set_d(x, P, GMP_RNDN);
    
    return *this;
}

my_mpfr& my_mpfr::operator = (const long double P)
{
    mpfr_set_ld(x, P, GMP_RNDN);
    
    return *this;
}

my_mpfr& my_mpfr::operator = (const int P)
{
    mpfr_set_si(x, (signed long int)P, GMP_RNDN);

    return *this;
}

my_mpfr& my_mpfr::operator = (const std::string a)
{
    mpfr_set_str(x, a.data(), 10, GMP_RNDN);

    return *this;
}

// Move constructor and move assignment

my_mpfr::my_mpfr(my_mpfr&& P)
{
    x->_mpfr_d = 0;
    //std::cout << "Calling move constructor my_mpfr(my_mpfr&&) for " << this << std::endl;
    mpfr_swap(x, P.x);
}

my_mpfr& my_mpfr::operator=(my_mpfr&& P)
{
    //std::cout << "Calling move assignment =(my_mpfr&&) for " << this << std::endl;
    mpfr_swap(x, P.x);
    return *this;
}

// Compund assignment operators

inline my_mpfr my_mpfr::operator+= (const my_mpfr& P)
{
    mpfr_add(x, x, P.x, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator+= (const unsigned long int P)
{
    mpfr_add_ui(x, x, P, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator+= (const signed long int P)
{
    mpfr_add_si(x, x, P, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator+= (const double P)
{
    mpfr_add_d(x, x, P, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator+= (const int P)
{
    mpfr_add_si(x, x,(signed long int)P, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator-= (const my_mpfr& P)
{
    mpfr_sub(x, x, P.x, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator-= (const unsigned long int P)
{
    mpfr_sub_ui(x, x, P, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator-= (const signed long int P)
{
    mpfr_sub_si(x, x, P, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator-= (const double P)
{
    mpfr_sub_d(x, x, P, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator-= (const int P)
{
    mpfr_sub_si(x, x,(signed long int)P, GMP_RNDN);

    return *this;
}


inline my_mpfr my_mpfr::operator*= (const my_mpfr& P)
{
    mpfr_mul(x, x, P.x, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator*= (const unsigned long int P)
{
    mpfr_mul_ui(x, x, P, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator*= (const signed long int P)
{
    mpfr_mul_si(x, x, P, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator*= (const double P)
{
    mpfr_mul_d(x, x, P, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator*= (const int P)
{
    mpfr_mul_si(x, x,(signed long int)P, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator/= (const my_mpfr& P)
{
    mpfr_div(x, x, P.x, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator/= (const unsigned long int P)
{
    mpfr_div_ui(x, x, P, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator/= (const signed long int P)
{
    mpfr_div_si(x, x, P, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator/= (const double P)
{
    mpfr_div_d(x, x, P, GMP_RNDN);

    return *this;
}

inline my_mpfr my_mpfr::operator/= (const int P)
{
    mpfr_div_si(x, x,(signed long int)P, GMP_RNDN);

    return *this;
}

// Arithmetic operators +,-,*,/

inline my_mpfr operator- (const my_mpfr& P)
{
    my_mpfr res;
    
    //mpfr_sub(res.x, my_mpfr(0.0).x, P.x, GMP_RNDN);
    
    mpfr_neg(res.x, P.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator+ (const my_mpfr& P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_add(res.x, P.x, Q.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator+ (const unsigned long int P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_add_ui(res.x, Q.x, P, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator+ (const my_mpfr& Q, const unsigned long int P)
{
    my_mpfr res;
    
    mpfr_add_ui(res.x, Q.x, P, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator+ (const signed long int P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_add_si(res.x, Q.x, P, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator+ (const my_mpfr& Q, const signed long int P)
{
    my_mpfr res;
    
    mpfr_add_si(res.x, Q.x, P, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator+ (const double P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_add_d(res.x, Q.x, P, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator+ (const my_mpfr& Q, const double P)
{
    my_mpfr res;
    
    mpfr_add_d(res.x, Q.x, P, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator+ (const int P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_add_si(res.x, Q.x, (signed long int)P, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator+ (const my_mpfr& Q, const int P)
{
    my_mpfr res;
    
    mpfr_add_si(res.x, Q.x, (signed long int)P, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator- (const my_mpfr& P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_sub(res.x, P.x, Q.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator- (const my_mpfr& P, const unsigned long int Q)
{
    my_mpfr res;
    
    mpfr_sub_ui(res.x, P.x, Q, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator- (const unsigned long int P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_ui_sub(res.x, P, Q.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator- (const my_mpfr& P, const signed long int Q)
{
    my_mpfr res;
    
    mpfr_sub_si(res.x, P.x, Q, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator- (const signed long int P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_si_sub(res.x, P, Q.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator- (const my_mpfr& P, const double Q)
{
    my_mpfr res;
    
    mpfr_sub_d(res.x, P.x, Q, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator- (const double P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_d_sub(res.x, P, Q.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator- (const my_mpfr& P, const int Q)
{
    my_mpfr res;
    
    mpfr_sub_si(res.x, P.x, (signed long int)Q, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator- (const int P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_si_sub(res.x, (signed long int)P, Q.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator* (const my_mpfr& P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_mul(res.x, P.x, Q.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator* (const my_mpfr& P, const unsigned long int Q)
{
    my_mpfr res;
    
    mpfr_mul_ui(res.x, P.x, Q, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator* (const unsigned long int P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_mul_ui(res.x, Q.x, P, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator* (const my_mpfr& P, const signed long int Q)
{
    my_mpfr res;
    
    mpfr_mul_si(res.x, P.x, Q, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator* (const signed long int P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_mul_si(res.x, Q.x, P, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator* (const my_mpfr& P, const double Q)
{
    my_mpfr res;
    
    mpfr_mul_d(res.x, P.x, Q, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator* (const double P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_mul_d(res.x, Q.x, P, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator* (const my_mpfr& P, const int Q)
{
    my_mpfr res;
    
    mpfr_mul_si(res.x, P.x, (signed long int)Q, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator* (const int P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_mul_si(res.x, Q.x, (signed long int)P, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator/ (const my_mpfr& P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_div(res.x, P.x, Q.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator/ (const my_mpfr& P, const unsigned long int Q)
{
    my_mpfr res;
    
    mpfr_div_ui(res.x, P.x, Q, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator/ (const unsigned long int P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_ui_div(res.x, P, Q.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator/ (const my_mpfr& P, const signed long int Q)
{
    my_mpfr res;
    
    mpfr_div_si(res.x, P.x, Q, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator/ (const signed long int P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_si_div(res.x, P, Q.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator/ (const my_mpfr& P, const double Q)
{
    my_mpfr res;
    
    mpfr_div_d(res.x, P.x, Q, GMP_RNDN);
    
    return res;
}
inline my_mpfr operator/ (const double P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_d_div(res.x, P, Q.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr operator/ (const my_mpfr& P, const int Q)
{
    my_mpfr res;
    
    mpfr_div_si(res.x, P.x, (signed long int)Q, GMP_RNDN);
    
    return res;
}
inline my_mpfr operator/ (const int P, const my_mpfr& Q)
{
    my_mpfr res;
    
    mpfr_si_div(res.x, (signed long int)P, Q.x, GMP_RNDN);
    
    return res;
}

// Relational operators

inline int operator== (const my_mpfr& P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp(P.x, Q.x);

    if(res == 0) res = 1;
    else res = 0;

    return res;
}

inline int operator!= (const my_mpfr& P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp(P.x, Q.x);

    if(res == 0) res = 0;
    else res = 1;

    return res;
}

inline int operator> (const my_mpfr& P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp(P.x, Q.x);

    if(res > 0) res = 1;
    else res = 0;

    return res;
}

inline int operator> (const my_mpfr& P, unsigned long int Q)
{
    int res;
    
    res = mpfr_cmp_ui(P.x, Q);
    
    if(res > 0) res = 1;
    else res = 0;

    return res;
}

inline int operator> (unsigned long int P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp_ui(Q.x, P);
    
    if(res < 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator> (const my_mpfr& P, signed long int Q)
{
    int res;
    
    res = mpfr_cmp_si(P.x, Q);
    
    if(res > 0) res = 1;
    else res = 0;

    return res;
}

inline int operator> (signed long int P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp_si(Q.x, P);
    
    if(res < 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator> (const my_mpfr& P, double Q)
{
    int res;
    
    res = mpfr_cmp_d(P.x, Q);
    
    if(res > 0) res = 1;
    else res = 0;

    return res;
}

inline int operator> (double P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp_d(Q.x, P);
    
    if(res < 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator< (const my_mpfr& P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp(P.x, Q.x);
    
    if(res < 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator< (const my_mpfr& P, unsigned long int Q)
{
    int res;
    
    res = mpfr_cmp_ui(P.x, Q);
    
    if(res < 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator< (unsigned long int P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp_ui(Q.x, P);
    
    if(res > 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator< (const my_mpfr& P, signed long int Q)
{
    int res;
    
    res = mpfr_cmp_si(P.x, Q);
    
    if(res < 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator< (signed long int P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp_si(Q.x, P);
    
    if(res > 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator< (const my_mpfr& P, double Q)
{
    int res;
    
    res = mpfr_cmp_d(P.x, Q);
    
    if(res < 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator< (double P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp_d(Q.x, P);
    
    if(res > 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator>= (const my_mpfr& P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp(P.x, Q.x);
    
    if(res >= 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator>= (const my_mpfr& P, unsigned long int Q)
{
    int res;
    
    res = mpfr_cmp_ui(P.x, Q);
    
    if(res >= 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator>= (unsigned long int P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp_ui(Q.x, P);
    
    if(res <= 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator>= (const my_mpfr& P, signed long int Q)
{
    int res;
    
    res = mpfr_cmp_si(P.x, Q);
    
    if(res >= 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator>= (signed long int P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp_si(Q.x, P);
    
    if(res <= 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator>= (const my_mpfr& P, double Q)
{
    int res;
    
    res = mpfr_cmp_d(P.x, Q);
    
    if(res >= 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator>= (double P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp_d(Q.x, P);
    
    if(res <= 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator<= (const my_mpfr& P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp(P.x, Q.x);
    
    if(res <= 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator<= (const my_mpfr& P, unsigned long int Q)
{
    int res;
    
    res = mpfr_cmp_ui(P.x, Q);
    
    if(res <= 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator<= (unsigned long int P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp_ui(Q.x, P);
    
    if(res >= 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator<= (const my_mpfr& P, signed long int Q)
{
    int res;
    
    res = mpfr_cmp_si(P.x, Q);
    
    if(res <= 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator<= (signed long int P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp_si(Q.x, P);
    
    if(res >= 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator<= (const my_mpfr& P, double Q)
{
    int res;
    
    res = mpfr_cmp_d(P.x, Q);
    
    if(res <= 0) res = 1;
    else res = 0;
    
    return res;
}

inline int operator<= (double P, const my_mpfr& Q)
{
    int res;
    
    res = mpfr_cmp_d(Q.x, P);
    
    if(res >= 0) res = 1;
    else res = 0;
    
    return res;
}

// Common functions

inline my_mpfr cos(const my_mpfr& P)
{
    my_mpfr res;
    
    mpfr_cos(res.x, P.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr sin(const my_mpfr& P)
{
    my_mpfr res;
    
    mpfr_sin(res.x, P.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr cosh(const my_mpfr& P)
{
    my_mpfr res;
    
    mpfr_cosh(res.x, P.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr sinh(const my_mpfr& P)
{
    my_mpfr res;
    
    mpfr_sinh(res.x, P.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr exp(const my_mpfr& P)
{
    my_mpfr res;
    
    mpfr_exp(res.x, P.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr exp2(const my_mpfr& P)
{
    my_mpfr res;
    
    mpfr_exp2(res.x, P.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr exp10(const my_mpfr& P)
{
    my_mpfr res;
    
    mpfr_exp10(res.x, P.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr log(const my_mpfr& P)
{
    my_mpfr res;
    
    mpfr_log(res.x, P.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr log2(const my_mpfr& P)
{
    my_mpfr res;
    
    mpfr_log2(res.x, P.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr log10(const my_mpfr& P)
{
    my_mpfr res;
    
    mpfr_log10(res.x, P.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr sqrt(const my_mpfr& P)
{
    my_mpfr res;
    
    mpfr_sqrt(res.x, P.x, GMP_RNDN);
 
    return res;
}

inline my_mpfr abs(const my_mpfr& P)
{
    my_mpfr res;
    
    mpfr_abs(res.x, P.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr pow(const my_mpfr& P, int N)
{
    my_mpfr res;
    
    res = my_mpfr(1.0);
    if(N>0){
	for(int i=0;i<N;i++) res = res*P;
    }
    else if (N<0){
	for(int i=0;i<-N;i++) res = res/P;
    }
    return res;
}

inline my_mpfr pow(const my_mpfr& P, const my_mpfr& x)
{
    my_mpfr res;
    
    mpfr_pow(res.x, P.x, x.x, GMP_RNDN);
    
    return res;
}

inline my_mpfr max(const my_mpfr& a, const my_mpfr& b)
{
    my_mpfr res;

    if (a>b) res=a;
    else res=b;

    return res;
}

inline my_mpfr max(const my_mpfr& a, const my_mpfr& b, const my_mpfr& c)
{
    my_mpfr res;

    res=max(a,b);
    if (c>res) res=c;

    return res;
}

inline my_mpfr max(const my_mpfr& a, const my_mpfr& b, const my_mpfr& c, const my_mpfr& d)
{
    my_mpfr res;

    res=max(a,b,c);
    if (d>res) res=d;

    return res;
}

inline my_mpfr round(const my_mpfr& P)
{
    my_mpfr res;

    mpfr_round(res.x, P.x);

    return res;
}

inline my_mpfr ceil(const my_mpfr& P)
{
    my_mpfr res;

    mpfr_ceil(res.x, P.x);

    return res;
}

inline my_mpfr floor(const my_mpfr& P) 
{
    my_mpfr res;

    mpfr_floor(res.x,P.x);

    return res;
}

inline int isfinite(const my_mpfr& P)
{
    return (!mpfr_nan_p(P.x) && !mpfr_inf_p(P.x));
}

inline int isinf(const my_mpfr& P)
{
    return mpfr_inf_p(P.x);
}

inline int isnan(const my_mpfr& P)
{
    return mpfr_nan_p(P.x);
}


std::ostream& operator<<(std::ostream& stream, const my_mpfr& op)
{
    static char str[100];
    char* str_ptr = str;
    mp_exp_t exp;
    mp_rnd_t rnd = mpfr_get_default_rounding_mode();
    mpfr_get_str(str, &exp, 10, 0, op.x, rnd);
    if ( str[0] == '-')
        stream << (str_ptr++)[0];
    stream << str_ptr[0]; 
    stream << '.';
    stream << (str_ptr+1);
    if (exp != 1) {
        stream << 'e';
        stream << exp-1;
    }
    return stream;
}


template<>
my_mpfr pi<my_mpfr>() {
    my_mpfr res;
    mpfr_const_pi(res.x,GMP_RNDN);
    return res;
}


#endif
