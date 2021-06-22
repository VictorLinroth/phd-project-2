
#ifndef INTERVAL_H
#define INTERVAL_H
/*

*/

#include<mpfi.h>
#include<mpfi_io.h>
#include"multiprec.hpp"

int MY_PREC = 256;

class INTERVAL
{
public:
	mpfi_t x;

    // Regular constructors and destructor
	INTERVAL();
    INTERVAL(const INTERVAL&);
	~INTERVAL();
    
	INTERVAL(const unsigned long int a);
	INTERVAL(const signed long int a);
	INTERVAL(const double a);
	INTERVAL(const int a);
    INTERVAL(const mpfr_t& a);
    INTERVAL(const my_mpfr& a);

	INTERVAL(const unsigned long int a, const unsigned long int b);
	INTERVAL(const signed long int a, const signed long int b);
	INTERVAL(double a, double b);
	INTERVAL(const int a, const int b);
    INTERVAL(const mpfr_t& a, const mpfr_t& b);
    INTERVAL(const my_mpfr& a, const my_mpfr& b);
    
    // Regular assignments
    INTERVAL& operator = (const INTERVAL&);
    INTERVAL& operator = (const unsigned long int);
    INTERVAL& operator = (const signed long int);
    INTERVAL& operator = (const double);
    INTERVAL& operator = (const int);
    INTERVAL& operator = (const mpfr_t&);
    INTERVAL& operator = (const my_mpfr&);

    // Move constructor and move assignment
    INTERVAL(INTERVAL&&);
    INTERVAL& operator=(INTERVAL&&);

	void setfull(double a, double b);
	void setfull(mpfr_t a, mpfr_t b);
	int position(double a);
	int position(mpfr_t a)const;
	int is_inside(double a);

	double left();
	double right();
	my_mpfr left_l() const;
	my_mpfr right_l() const;
    INTERVAL mid();

    // Compound assigment operators
    inline INTERVAL operator+= (const INTERVAL&);
    inline INTERVAL operator+= (const unsigned long int);
    inline INTERVAL operator+= (const signed long int);
    inline INTERVAL operator+= (const double);
    inline INTERVAL operator+= (const int);
    inline INTERVAL operator+= (const my_mpfr&);

    inline INTERVAL operator-= (const INTERVAL&);
    inline INTERVAL operator-= (const unsigned long int);
    inline INTERVAL operator-= (const signed long int);
    inline INTERVAL operator-= (const double);
    inline INTERVAL operator-= (const int);
    inline INTERVAL operator-= (const my_mpfr&);

    inline INTERVAL operator*= (const INTERVAL&);
    inline INTERVAL operator*= (const unsigned long int);
    inline INTERVAL operator*= (const signed long int);
    inline INTERVAL operator*= (const double);
    inline INTERVAL operator*= (const int);
    inline INTERVAL operator*= (const my_mpfr&);

    inline INTERVAL operator/= (const INTERVAL&);
    inline INTERVAL operator/= (const unsigned long int);
    inline INTERVAL operator/= (const signed long int);
    inline INTERVAL operator/= (const double);
    inline INTERVAL operator/= (const int);
    inline INTERVAL operator/= (const my_mpfr&);

    // Efficient swap function
    friend void swap(INTERVAL& P, INTERVAL& Q)
    {
        mpfi_swap(P.x, Q.x);
    }
};

// Arithemtic operations +,-,*,/
inline INTERVAL operator- (const INTERVAL& P);

inline INTERVAL operator+ (const INTERVAL&, const INTERVAL&);
inline INTERVAL operator+ (const unsigned long int, const INTERVAL&);
inline INTERVAL operator+ (const INTERVAL&, const unsigned long int);
inline INTERVAL operator+ (const signed long int, const INTERVAL&);
inline INTERVAL operator+ (const INTERVAL&, const signed long int);
inline INTERVAL operator+ (const double, const INTERVAL&);
inline INTERVAL operator+ (const INTERVAL&, const double);
inline INTERVAL operator+ (const int, const INTERVAL&);
inline INTERVAL operator+ (const INTERVAL&, const int);
inline INTERVAL operator+ (const mpfr_t&, const INTERVAL&);
inline INTERVAL operator+ (const INTERVAL&, const mpfr_t&);
inline INTERVAL operator+ (const my_mpfr&, const INTERVAL&);
inline INTERVAL operator+ (const INTERVAL&, const my_mpfr&);

inline INTERVAL operator- (const INTERVAL&, const INTERVAL&);
inline INTERVAL operator- (const unsigned long int, const INTERVAL&);
inline INTERVAL operator- (const INTERVAL&, const unsigned long int);
inline INTERVAL operator- (const signed long int, const INTERVAL&);
inline INTERVAL operator- (const INTERVAL&, const signed long int);
inline INTERVAL operator- (const double, const INTERVAL&);
inline INTERVAL operator- (const INTERVAL&, const double);
inline INTERVAL operator- (const int, const INTERVAL&);
inline INTERVAL operator- (const INTERVAL&, const int);
inline INTERVAL operator- (const mpfr_t&, const INTERVAL&);
inline INTERVAL operator- (const INTERVAL&, const mpfr_t&);
inline INTERVAL operator- (const my_mpfr&, const INTERVAL&);
inline INTERVAL operator- (const INTERVAL&, const my_mpfr&);

inline INTERVAL operator* (const INTERVAL&, const INTERVAL&);
inline INTERVAL operator* (const unsigned long int, const INTERVAL&);
inline INTERVAL operator* (const INTERVAL&, const unsigned long int);
inline INTERVAL operator* (const signed long int, const INTERVAL&);
inline INTERVAL operator* (const INTERVAL&, const signed long int);
inline INTERVAL operator* (const double, const INTERVAL&);
inline INTERVAL operator* (const INTERVAL&, const double);
inline INTERVAL operator* (const int, const INTERVAL&);
inline INTERVAL operator* (const INTERVAL&, const int);
inline INTERVAL operator* (const mpfr_t&, const INTERVAL&);
inline INTERVAL operator* (const INTERVAL&, const mpfr_t&);
inline INTERVAL operator* (const my_mpfr&, const INTERVAL&);
inline INTERVAL operator* (const INTERVAL&, const my_mpfr&);

inline INTERVAL operator/ (const INTERVAL&, const INTERVAL&);
inline INTERVAL operator/ (const unsigned long int, const INTERVAL&);
inline INTERVAL operator/ (const INTERVAL&, const unsigned long int);
inline INTERVAL operator/ (const signed long int, const INTERVAL&);
inline INTERVAL operator/ (const INTERVAL&, const signed long int);
inline INTERVAL operator/ (const double, const INTERVAL&);
inline INTERVAL operator/ (const INTERVAL&, const double);
inline INTERVAL operator/ (const int, const INTERVAL&);
inline INTERVAL operator/ (const INTERVAL&, const int);
inline INTERVAL operator/ (const mpfr_t&, const INTERVAL&);
inline INTERVAL operator/ (const INTERVAL&, const mpfr_t&);
inline INTERVAL operator/ (const my_mpfr&, const INTERVAL&);
inline INTERVAL operator/ (const INTERVAL&, const my_mpfr&);

inline INTERVAL abs(const INTERVAL& P);
inline INTERVAL sin(const INTERVAL& P);
inline INTERVAL sinh(const INTERVAL& P);
inline INTERVAL cos(const INTERVAL& P);
inline INTERVAL cosh(const INTERVAL& P);
inline INTERVAL exp(const INTERVAL& P);
inline INTERVAL exp2(const INTERVAL& P);
inline INTERVAL log(const INTERVAL& P);
inline INTERVAL log2(const INTERVAL& P);
inline INTERVAL sqr(const INTERVAL& P);
inline INTERVAL sqrt(const INTERVAL& P);
inline INTERVAL pow(const INTERVAL& P, int N);
inline INTERVAL pow(const INTERVAL& P, const INTERVAL& Q);
inline INTERVAL max(const INTERVAL& P, const INTERVAL& Q);

inline INTERVAL const_pi();

inline INTERVAL convhull(INTERVAL, const INTERVAL&);

INTERVAL str_to_INTERVAL(char *s);

bool check_invertibility(INTERVAL **A, int dim);

// Stream operator overload
std::ostream& operator<<(std::ostream& stream, const INTERVAL& op);

/*-------------------------------------------------------------------------------------------------------------------------------------------*/

INTERVAL::INTERVAL () 
{
	mpfi_init_set_d(x, 0.);
}

INTERVAL::INTERVAL(const INTERVAL& a)
{
    mpfi_init_set(x, a.x);
}

INTERVAL::INTERVAL(const unsigned long int a)
{
    mpfi_init_set_ui(x,a);
}

INTERVAL::INTERVAL(const signed long int a)
{
    mpfi_init_set_si(x,a);
}

INTERVAL::INTERVAL (double a) 
{
	mpfi_init_set_d(x, a);
}

INTERVAL::INTERVAL (int a) 
{
	mpfi_init_set_si(x, (int)a);
}

INTERVAL::INTERVAL (const mpfr_t& a) 
{
	mpfi_init_set_fr(x, a);
}

INTERVAL::INTERVAL (const my_mpfr& a) 
{
	mpfi_init_set_fr(x, a.x);
}

INTERVAL::INTERVAL(const unsigned long int a, const unsigned long int b)
{
    mpfi_init(x);
    mpfi_interv_ui(x, a, b);
}

INTERVAL::INTERVAL(const signed long int a, const signed long int b)
{
    mpfi_init(x);
    mpfi_interv_si(x, a, b);
}

INTERVAL::INTERVAL (double a, double b) 
{
	mpfi_init(x);
	mpfi_interv_d(x, a, b);
}

INTERVAL::INTERVAL(const int a, const int b)
{
    mpfi_init(x);
    mpfi_interv_si(x, (int)a, (int)b);
}

INTERVAL::INTERVAL(const mpfr_t& a, const mpfr_t& b)
{
    mpfi_init(x);
    mpfi_interv_fr(x, a, b);
}

INTERVAL::INTERVAL(const my_mpfr& a, const my_mpfr& b)
{
    mpfi_init(x);
    mpfi_interv_fr(x, a.x, b.x);
}


INTERVAL::INTERVAL(INTERVAL&& P)
{
    x->left._mpfr_d = 0;
    mpfi_swap(x, P.x);
}

INTERVAL& INTERVAL::operator=(INTERVAL&& P)
{
    mpfi_swap(x, P.x);
    return *this;
}

INTERVAL::~INTERVAL()
{
//    std::cout << "Destructor for INTERVAL = " << *this << std::endl;
    if (x->left._mpfr_d!=0)
        mpfi_clear(x);
}

INTERVAL& INTERVAL::operator = (const INTERVAL& P)
{
	if(this != &P)
	{
		mpfi_set(x, P.x);
	}
	return *this;
}

INTERVAL& INTERVAL::operator = (const unsigned long int P)
{
	mpfi_set_ui(x, P);
	return *this;
}

INTERVAL& INTERVAL::operator = (const signed long int P)
{
	mpfi_set_si(x, P);
	return *this;
}

INTERVAL& INTERVAL::operator = (const double P)
{
	mpfi_set_d(x, P);
	return *this;
}

INTERVAL& INTERVAL::operator = (const int P)
{
	mpfi_set_si(x, (int)P);
	return *this;
}

INTERVAL& INTERVAL::operator = (const mpfr_t& P)
{
	mpfi_set_fr(x, P);
	return *this;
}

INTERVAL& INTERVAL::operator = (const my_mpfr& P)
{
	mpfi_set_fr(x, P.x);
	return *this;
}

void INTERVAL::setfull(double a, double b)
{
	mpfi_interv_d(x, a, b);

	return ;
}

void INTERVAL::setfull(mpfr_t a, mpfr_t b)
{
	mpfi_interv_fr(x, a, b);

	return ;
}

int INTERVAL::position(double a)
{
	mpfr_t left, right;
	int flag = 0, retl, retr;

	mpfr_init(left);
	mpfr_init(right);

	mpfi_get_left(left, x);
	mpfi_get_right(right, x);

	retl = mpfr_cmp_d(left, a);
	retr = mpfr_cmp_d(right, a);
	
	if(retl > 0)
	{
		flag = -1;
	}
	if(retr < 0)
	{
		flag = 1;
	}

	mpfr_clear(left);
	mpfr_clear(right);

	return flag;
}

int INTERVAL::position(mpfr_t a)const
{
	mpfr_t left, right;
	int flag = 0, retl, retr;

	mpfr_init(left);
	mpfr_init(right);

	mpfi_get_left(left, x);
	mpfi_get_right(right, x);

	retl = mpfr_cmp(left, a);
	retr = mpfr_cmp(right, a);
	
	if(retl > 0)
	{
		flag = -1;
	}
	if(retr < 0)
	{
		flag = 1;
	}

	mpfr_clear(left);
	mpfr_clear(right);

	return flag;
}

//Returns 0 if a is not inside x, 1 otherwise
int INTERVAL::is_inside(double a)
{
	int flag;

	flag = mpfi_is_inside_d(a, x);

	return flag;
}

double INTERVAL::left(void )
{	
	mpfr_t left;
	double res;

	mpfr_init(left);

	mpfi_get_left(left, x);
	res = mpfr_get_d(left, GMP_RNDD);

	mpfr_clear(left);

	return res;
}

double INTERVAL::right(void )
{	
	mpfr_t right;
	double res;

	mpfr_init(right);

	mpfi_get_right(right, x);
	res = mpfr_get_d(right, GMP_RNDU);

	mpfr_clear(right);

	return res;
}

my_mpfr INTERVAL::left_l(void ) const
{	
   my_mpfr res;

	mpfi_get_left(res.x, x);

	return res;
}

my_mpfr INTERVAL::right_l(void ) const
{	
   my_mpfr res;

	mpfi_get_right(res.x, x);

	return res;
}

INTERVAL INTERVAL::mid(void)
{	
	mpfr_t left;
   mpfr_t right;
   INTERVAL l, r, res;

	mpfr_init(left);
	mpfr_init(right);

	mpfi_get_left(left, x);
	mpfi_get_right(right, x);
   
   mpfi_set_fr(l.x, left);
   mpfi_set_fr(r.x, right);
  
   res = (l+r)/INTERVAL(2.);

	mpfr_clear(left);
	mpfr_clear(right);

	return res;
}

// Compound assigment operators
inline INTERVAL INTERVAL::operator+= (const INTERVAL& P)
{
    mpfi_add(x, x, P.x);

    return *this;
}

inline INTERVAL INTERVAL::operator+= (const unsigned long int P)
{
    mpfi_add_ui(x, x, P);

    return *this;
}

inline INTERVAL INTERVAL::operator+= (const signed long int P)
{
    mpfi_add_si(x, x, P);

    return *this;
}

inline INTERVAL INTERVAL::operator+= (const double P)
{
    mpfi_add_d(x, x, P);

    return *this;
}

inline INTERVAL INTERVAL::operator+= (const int P)
{
    mpfi_add_si(x, x, (signed long int)P);

    return *this;
}

inline INTERVAL INTERVAL::operator+= (const my_mpfr& P)
{
    mpfi_add_fr(x, x, P.x);

    return *this;
}

inline INTERVAL INTERVAL::operator-= (const INTERVAL& P)
{
    mpfi_sub(x, x, P.x);

    return *this;
}

inline INTERVAL INTERVAL::operator-= (const unsigned long int P)
{
    mpfi_sub_ui(x, x, P);

    return *this;
}

inline INTERVAL INTERVAL::operator-= (const signed long int P)
{
    mpfi_sub_si(x, x, P);

    return *this;
}

inline INTERVAL INTERVAL::operator-= (const double P)
{
    mpfi_sub_d(x, x, P);

    return *this;
}

inline INTERVAL INTERVAL::operator-= (const int P)
{
    mpfi_sub_si(x, x, (signed long int)P);

    return *this;
}

inline INTERVAL INTERVAL::operator-= (const my_mpfr& P)
{
    mpfi_sub_fr(x, x, P.x);

    return *this;
}

inline INTERVAL INTERVAL::operator*= (const INTERVAL& P)
{
    mpfi_mul(x, x, P.x);

    return *this;
}

inline INTERVAL INTERVAL::operator*= (const unsigned long int P)
{
    mpfi_mul_ui(x, x, P);

    return *this;
}

inline INTERVAL INTERVAL::operator*= (const signed long int P)
{
    mpfi_mul_si(x, x, P);

    return *this;
}

inline INTERVAL INTERVAL::operator*= (const double P)
{
    mpfi_mul_d(x, x, P);

    return *this;
}

inline INTERVAL INTERVAL::operator*= (const int P)
{
    mpfi_mul_si(x, x, (signed long int)P);

    return *this;
}

inline INTERVAL INTERVAL::operator*= (const my_mpfr& P)
{
    mpfi_mul_fr(x, x, P.x);

    return *this;
}

inline INTERVAL INTERVAL::operator/= (const INTERVAL& P)
{
    mpfi_div(x, x, P.x);

    return *this;
}

inline INTERVAL INTERVAL::operator/= (const unsigned long int P)
{
    mpfi_div_ui(x, x, P);

    return *this;
}

inline INTERVAL INTERVAL::operator/= (const signed long int P)
{
    mpfi_div_si(x, x, P);

    return *this;
}

inline INTERVAL INTERVAL::operator/= (const double P)
{
    mpfi_div_d(x, x, P);

    return *this;
}

inline INTERVAL INTERVAL::operator/= (const int P)
{
    mpfi_div_si(x, x, (signed long int)P);

    return *this;
}

inline INTERVAL INTERVAL::operator/= (const my_mpfr& P)
{
    mpfi_div_fr(x, x, P.x);

    return *this;
}


// Arithemtic operations +,-,*,/

inline INTERVAL operator- (const INTERVAL& P)
{
	INTERVAL res;

	//mpfi_add(res.x, INTERVAL(0).x, P.x);

    mpfi_neg(res.x, P.x);

	return res;
}

inline INTERVAL operator+ (const INTERVAL& P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_add(res.x, P.x, Q.x);

	return res;
}

inline INTERVAL operator+ (const unsigned long int P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_add_ui(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator+ (const INTERVAL& Q, const unsigned long int P)
{
	INTERVAL res;

	mpfi_add_ui(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator+ (const signed long int P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_add_si(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator+ (const INTERVAL& Q, const signed long int P)
{
	INTERVAL res;

	mpfi_add_si(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator+ (const double P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_add_d(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator+ (const INTERVAL& Q, const double P)
{
	INTERVAL res;

	mpfi_add_d(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator+ (const int P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_add_si(res.x, Q.x, (signed long int)P);

	return res;
}

inline INTERVAL operator+ (const INTERVAL& Q, const int P)
{
	INTERVAL res;

	mpfi_add_si(res.x, Q.x, (signed long int)P);

	return res;
}

inline INTERVAL operator+ (const mpfr_t& P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_add_fr(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator+ (const INTERVAL& Q, const mpfr_t& P)
{
	INTERVAL res;

	mpfi_add_fr(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator+ (const my_mpfr& P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_add_fr(res.x, Q.x, P.x);

	return res;
}

inline INTERVAL operator+ (const INTERVAL& Q, const my_mpfr& P)
{
	INTERVAL res;

	mpfi_add_fr(res.x, Q.x, P.x);

	return res;
}

inline INTERVAL operator- (const INTERVAL& P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_sub(res.x, P.x, Q.x);

	return res;
}

inline INTERVAL operator- (const INTERVAL& P, const unsigned long int Q)
{
	INTERVAL res;

	mpfi_sub_ui(res.x, P.x, Q);

	return res;
}

inline INTERVAL operator- (const unsigned long int P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_ui_sub(res.x, P, Q.x);

	return res;
}

inline INTERVAL operator- (const INTERVAL& P, const signed long int Q)
{
	INTERVAL res;

	mpfi_sub_si(res.x, P.x, Q);

	return res;
}

inline INTERVAL operator- (const signed long int P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_si_sub(res.x, P, Q.x);

	return res;
}

inline INTERVAL operator- (const INTERVAL& P, const double Q)
{
	INTERVAL res;

	mpfi_sub_d(res.x, P.x, Q);

	return res;
}

inline INTERVAL operator- (const double P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_d_sub(res.x, P, Q.x);

	return res;
}

inline INTERVAL operator- (const INTERVAL& P, const int Q)
{
	INTERVAL res;

	mpfi_sub_si(res.x, P.x, Q);

	return res;
}

inline INTERVAL operator- (const int P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_si_sub(res.x, P, Q.x);

	return res;
}

inline INTERVAL operator- (const INTERVAL& P, const mpfr_t& Q)
{
	INTERVAL res;

	mpfi_sub_fr(res.x, P.x, Q);

	return res;
}

inline INTERVAL operator- (const mpfr_t& P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_fr_sub(res.x, P, Q.x);

	return res;
}

inline INTERVAL operator- (const INTERVAL& P, const my_mpfr& Q)
{
	INTERVAL res;

	mpfi_sub_fr(res.x, P.x, Q.x);

	return res;
}

inline INTERVAL operator- (const my_mpfr& P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_fr_sub(res.x, P.x, Q.x);

	return res;
}

inline INTERVAL operator* (const INTERVAL& P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_mul(res.x, P.x, Q.x);

	return res;
}

inline INTERVAL operator* (const unsigned long int P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_mul_ui(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator* (const INTERVAL& Q, const unsigned long int P)
{
	INTERVAL res;

	mpfi_mul_ui(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator* (const signed long int P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_mul_si(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator* (const INTERVAL& Q, const signed long int P)
{
	INTERVAL res;

	mpfi_mul_si(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator* (const double P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_mul_d(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator* (const INTERVAL& Q, const double P)
{
	INTERVAL res;

	mpfi_mul_d(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator* (const int P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_mul_si(res.x, Q.x, (signed long int)P);

	return res;
}

inline INTERVAL operator* (const INTERVAL& Q, const int P)
{
	INTERVAL res;

	mpfi_mul_si(res.x, Q.x, (signed long int)P);

	return res;
}

inline INTERVAL operator* (const mpfr_t& P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_mul_fr(res.x, Q.x, P);

	return res;
}

inline INTERVAL operator* (const INTERVAL& P, const mpfr_t& Q)
{
	INTERVAL res;

	mpfi_mul_fr(res.x, P.x, Q);

	return res;
}

inline INTERVAL operator* (const my_mpfr& P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_mul_fr(res.x, Q.x, P.x);

	return res;
}

inline INTERVAL operator* (const INTERVAL& P, const my_mpfr& Q)
{
	INTERVAL res;

	mpfi_mul_fr(res.x, P.x, Q.x);

	return res;
}

inline INTERVAL operator/ (const INTERVAL& P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_div(res.x, P.x, Q.x);

	return res;
}

inline INTERVAL operator/ (const INTERVAL& P, const unsigned long int Q)
{
	INTERVAL res;

	mpfi_div_ui(res.x, P.x, Q);

	return res;
}

inline INTERVAL operator/ (const unsigned long int P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_ui_div(res.x, P, Q.x);

	return res;
}

inline INTERVAL operator/ (const INTERVAL& P, const signed long int Q)
{
	INTERVAL res;

	mpfi_div_si(res.x, P.x, Q);

	return res;
}

inline INTERVAL operator/ (const signed long int P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_si_div(res.x, P, Q.x);

	return res;
}

inline INTERVAL operator/ (const INTERVAL& P, const double Q)
{
	INTERVAL res;

	mpfi_div_d(res.x, P.x, Q);

	return res;
}

inline INTERVAL operator/ (const double P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_d_div(res.x, P, Q.x);

	return res;
}

inline INTERVAL operator/ (const INTERVAL& P, const int Q)
{
	INTERVAL res;

	mpfi_div_si(res.x, P.x, Q);

	return res;
}

inline INTERVAL operator/ (const int P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_si_div(res.x, P, Q.x);

	return res;
}

inline INTERVAL operator/ (const INTERVAL& P, const mpfr_t& Q)
{
	INTERVAL res;

	mpfi_div_fr(res.x, P.x, Q);

	return res;
}

inline INTERVAL operator/ (const mpfr_t& P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_fr_div(res.x, P, Q.x);

	return res;
}

inline INTERVAL operator/ (const INTERVAL& P, const my_mpfr& Q)
{
	INTERVAL res;

	mpfi_div_fr(res.x, P.x, Q.x);

	return res;
}

inline INTERVAL operator/ (const my_mpfr& P, const INTERVAL& Q)
{
	INTERVAL res;

	mpfi_fr_div(res.x, P.x, Q.x);

	return res;
}

inline INTERVAL abs(const INTERVAL& P)
{
	INTERVAL res;

	mpfi_abs(res.x, P.x);

	return res;
}

inline INTERVAL cos(const INTERVAL& P)
{
	INTERVAL res;

	mpfi_cos(res.x, P.x);

	return res;
}

inline INTERVAL sin(const INTERVAL& P)
{
	INTERVAL res;

	mpfi_sin(res.x, P.x);

	return res;
}

inline INTERVAL cosh(const INTERVAL& P)
{
	INTERVAL res;

	mpfi_cosh(res.x, P.x);

	return res;
}

inline INTERVAL sinh(const INTERVAL& P)
{
	INTERVAL res;

	mpfi_sinh(res.x, P.x);

	return res;
}

inline INTERVAL exp(const INTERVAL& P)
{
	INTERVAL res;

	mpfi_exp(res.x, P.x);

	return res;
}

inline INTERVAL exp2(const INTERVAL& P)
{
	INTERVAL res;

	mpfi_exp2(res.x, P.x);

	return res;
}

inline INTERVAL log(const INTERVAL& P)
{
	INTERVAL res;

	mpfi_log(res.x, P.x);

	return res;
}

inline INTERVAL log2(const INTERVAL& P)
{
	INTERVAL res;

	mpfi_log2(res.x, P.x);

	return res;
}

inline INTERVAL sqr(const INTERVAL& P)
{
	INTERVAL res;

	mpfi_sqr(res.x, P.x);

	return res;
}

inline INTERVAL sqrt(const INTERVAL& P)
{
	INTERVAL res;

	mpfi_sqrt(res.x, P.x);

	return res;
}

inline INTERVAL pow(const INTERVAL& P, int N)
{
	INTERVAL res;
	int i;

	res = INTERVAL(1);
   if(N > 0)
   {
      for(i = 0; i < N; i++)
      {
         res = res*P;
      }
   }
   if(N < 0)
   {
      for(i = 0; i < -N; i++)
      {
         res = res/P;
      }
   }

	return res;
}

inline INTERVAL pow(const INTERVAL& P, const INTERVAL& Q)
{
    return exp(Q*log(P));
}

inline INTERVAL max(const INTERVAL& P, const INTERVAL& Q)
{
    INTERVAL res(max(P.left_l(),Q.left_l()),max(P.right_l(),Q.right_l()));
    return res;
}

inline INTERVAL convhull(INTERVAL P, const INTERVAL& Q)
{
    mpfi_put(P.x, Q.x);

    return P;
}

inline INTERVAL const_pi()
{
   INTERVAL res;

   mpfi_const_pi(res.x);

   return res;
}

INTERVAL str_to_INTERVAL(char *s)
{
	INTERVAL res;

	mpfi_set_str(res.x, s, 10);

	return res;
}

bool check_invertibility(INTERVAL **A, int dim)
{
	bool res = true;
	int i, j, k, flag = 0, ipivotpos = 0, jpivotpos = 0;
	INTERVAL pivot, pivotline, pivotabs, *aux, aux0, aux1;
	mpfr_t left0, left1;

	mpfr_init(left0);
	mpfr_init(left1);

	aux = new INTERVAL[dim];

	for(i = 0; i < dim; i++)
	{
	/*	for(j = 0; j < dim; j++)
		{
			for(k = 0; k < dim; k++)
			{
				mpfi_out_str(stdout, 10, 8, A[j][k].x);
			}
			printf("\n");
		}
		printf("\n");*/
		flag = 0;
		for(j = i; j < dim; j++)
		{
			for(k = j; k < dim; k++)
			{
				mpfi_abs(aux0.x, A[j][k].x);
				if(mpfi_has_zero(aux0.x) > 0)
				{
					continue;
				}else
				{
					if(flag == 0)
					{
						ipivotpos = j;
						jpivotpos = k;
						pivotabs = aux0;
						flag = 1;
					}else
					{
						mpfi_mig(left0, aux0.x);
						mpfi_mig(left1, pivotabs.x);
						if(mpfr_cmp(left0, left1) > 0)
						{
							pivotabs = aux0;
							ipivotpos = j;
							jpivotpos = k;
						}
					}
				}
			}
		}
		if(flag == 0)
		{
			res = false;
			break;
		}
		if(ipivotpos != i)
		{
			for(j = 0; j < dim; j++)
			{
				aux[j] = A[ipivotpos][j];
				A[ipivotpos][j] = A[i][j];
				A[i][j] = aux[j];
			}
		}
		if(jpivotpos != i)
		{
			for(j = 0; j < dim; j++)
			{
				aux[j] = A[j][jpivotpos];
				A[j][jpivotpos] = A[j][i];
				A[j][i] = aux[j];
			}
		}
		pivot = A[i][i];
		for(j = i+1; j < dim; j++)
		{
			pivotline = A[j][i]/pivot;
			for(k = i; k < dim; k++)
			{
				A[j][k] = A[j][k]-A[i][k]*pivotline;
			}
		}
	}

	delete []aux;

	return res;
}

/*inline INTERVAL abs(const INTERVAL& P)
{
        INTERVAL res;
        mpfr_t val0, leftP, rightP, aux;

        mpfr_init(val0);
        mpfr_init(leftP);
        mpfr_init(rightP);
        mpfr_init(aux);

        mpfr_set_d(val0, 0.0, GMP_RNDN);
        mpfi_get_left(leftP, P.x);
        mpfi_get_right(rightP, P.x);
        
        int a = P.position(val0);
        if(P.position(val0)==-1){
            res=P;
        }
        else if (P.position(val0)==0){
            mpfr_sub(aux,val0,leftP, GMP_RNDN);
            if (mpfr_cmp(rightP, aux)){
                res.setfull(val0,rightP);
            }
            else{
                res.setfull(val0,aux);
            }
        }
        else if (P.position(val0)==1){
            mpfr_sub(aux,val0,leftP, GMP_RNDN);
            mpfr_sub(leftP,val0,rightP, GMP_RNDN);
            res.setfull(leftP,aux);
        }

        mpfr_clear(val0);
        mpfr_clear(leftP);
        mpfr_clear(rightP);
        mpfr_clear(aux);

        return res;
}
*/

std::ostream& operator<<(std::ostream& stream, const INTERVAL& op) {
        stream << "[" << op.left_l() << "," << op.right_l() << "]";
        return stream;
}

template<>
INTERVAL pi<INTERVAL>() {
    INTERVAL res;
    mpfi_const_pi(res.x);
    return res;
}



#endif
