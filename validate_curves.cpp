#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <algorithm>
#include "Interval.hpp"
#include "TransformTools.hpp"
#include "DSTools.hpp"

const std::size_t dim = 2;
const std::size_t prec = 200;


using Real = my_mpfr;
using Vec2 = Vec<my_mpfr,dim>;
using Vec2_int = Vec<INTERVAL,dim>;

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
    printf("Using my_mpfr with precision %lu (MPFR_PREC_MIN = %d, MPFR_PREC_MAX = %ld)\n",prec,MPFR_PREC_MIN,MPFR_PREC_MAX);
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
    //const ValueType twopirec = ValueType(1)/twopi;
    x = param.lambda*y + param.mu + param.epsilon/twopi*sin(twopi*x);
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

INTERVAL hurwitz_zeta_function(int s, INTERVAL q, int terms) {
    INTERVAL res(0);

    for (int n=0; n!=terms; ++n)
        res += pow(n+q,-s);
    res += pow(terms+q,1-s)/(s-1);

    return res;
}

int gamma_bound(INTERVAL x) {
    int res(1);
    for(int n=1; INTERVAL(n).right_l() < x.right_l(); ++n)
        res *= n;
    return res;
}

int main(void) {

    init<Real>();
    
    const INTERVAL twopi = 2*pi<INTERVAL>();

    // -------- Calculate interval for omega --------------
    INTERVAL omega = (sqrt(INTERVAL(5.0))-INTERVAL(1.0))/INTERVAL(2.0);

    Real omega_real = omega.left_l() + (omega.right_l() - omega.left_l())/2;

    std::cout << "omega = " << omega << std::endl;
    
    std::vector<std::size_t> NumCoeff_list = {1024, 2048, 4096};
    //std::size_t NumCoeff = NumCoeff_list[0];
    std::size_t NumCoeff_approx =128; //240;
    std::size_t OrbitSize = (1<<15) - NumCoeff_approx;
    std::size_t ConvIter = 10000; //

    //std::cout << "NumCoeff= " << NumCoeff<< std::endl;
    std::cout << "NumCoeff_approx = " << NumCoeff_approx << std::endl;
    std::cout << "Using OrbitSize: " << OrbitSize << std::endl;
    std::cout << "Using ConvIter: " << ConvIter << std::endl;

    std::vector<Vec2> K_p_approx(2*NumCoeff_approx);
    std::vector<Vec2_int> K_p_Coeff;
    std::vector<Vec2_int> K_p_ug;
    std::vector<Vec2_int> E_Coeff;
    std::vector<Vec2_int> DK_Coeff;
    std::vector<Vec2_int> DK_omega;
    std::vector<Vec2_int> N0;
    std::vector<Vec2_int> N0_tilde;
    std::vector<Vec2_int> N0_omega;
    std::vector<INTERVAL> htilde;
    std::vector<INTERVAL> Ttilde;
    std::vector<INTERVAL> BLtilde;
    std::vector<INTERVAL> BN;
    std::vector<INTERVAL> RBN;
    std::vector<INTERVAL> tempvector;
    std::vector<INTERVAL> Im_K;

    std::vector<INTERVAL> rho_list(10);
    rho_list[0] = my_mpfr(1)/1024;
    for( std::size_t i = 1; i!=rho_list.size(); ++i)
        rho_list[i] = rho_list[i-1].right_l()/2;
    INTERVAL rho = 0.001;
    INTERVAL rho_hat = 0.01;
    INTERVAL rho_infinity = rho/4;
    INTERVAL delta = rho/4;


    INTERVAL gamma = 1.38, tau = 1.3;

    std::cout << "gamma = " << gamma.right_l() << "\ntau = " << tau.right_l() << std::endl;



    const int d = 1;

    INTERVAL c_R = sqrt(pow(INTERVAL(2),d-3)*hurwitz_zeta_function(2,pow(INTERVAL(2),tau),100)*gamma_bound(2*tau+1))/pow(twopi,tau); // Constant for Russmann esitmate

    std::cout << "c_R = " << c_R.right_l() << std::endl;

    std::fstream parameter_list("parameters.txt");
    //std::fstream parameter_list("parameter_list.txt");

    if (!parameter_list.is_open()) {
        std::cout << "Could not open file!" << std::endl;
        return 0;
    }

    while (!parameter_list.eof()) {


        // ------- Read parameters from trace_curve -----------
        std::string mu_str, epsilon_str, lambda_str;
        parameter_list >> mu_str;
        parameter_list >> epsilon_str;
        parameter_list >> lambda_str;
        MapParameters<INTERVAL> MP_int = { INTERVAL(mu_str), INTERVAL(epsilon_str), INTERVAL(lambda_str)};
        MapParameters<Real> MP_real = { Real(mu_str), Real(epsilon_str), Real(lambda_str)};
        std::cout << "mu = " << MP_real.mu << ", epsilon = " << MP_real.epsilon << ", lambda = " << MP_real.lambda << std::endl;
        
        INTERVAL lambda = MP_int.lambda, mu = MP_int.mu, epsilon = MP_int.epsilon;

        // -------- Make Orbit ---------------------
   
        { // Create scope for K_p_Orbit
            std::vector<Vec2> K_p_Orbit(2*OrbitSize);

            Vec2 z = {Real(0),omega_real};
            //z[1] = omega_real;

            for( std::size_t i=0; i!=ConvIter; ++i)
                z = F(z,MP_real);

            Real zphase = (OrbitSize/2)*omega_real;
            zphase -= floor(zphase);

            K_p_Orbit[0][0] = z[0] + zphase;
            K_p_Orbit[0][1] = z[1];
            for( int i=-OrbitSize/2+1; i!=OrbitSize/2; ++i) {
                K_p_Orbit[2*(i+OrbitSize/2)][1] = phi(z[0],z[1],MP_real);
                K_p_Orbit[2*(i+OrbitSize/2)][0] = K_p_Orbit[2*(i+OrbitSize/2)][1];
                K_p_Orbit[2*(i+OrbitSize/2)][0] += K_p_Orbit[2*(i-1+OrbitSize/2)][0] - omega_real;
                z[0] = i*omega_real + K_p_Orbit[2*(i+OrbitSize/2)][0];
                z[0] -= floor(z[0]);
                z[1] = K_p_Orbit[2*(i+OrbitSize/2)][1];
            }

            // -------- Calculate Fourier Coefficients ------------

            fourier_coeff_qg(K_p_approx, K_p_Orbit, omega_real);

            K_p_approx[NumCoeff_approx] = K_p_approx[NumCoeff_approx+1] = Vec2(0); // Since K_p_Orbit is real
        }

        bool existence = false;

        for( std::size_t n=0; n!=NumCoeff_list.size() && !existence; ++n) {

            std::size_t NumCoeff = NumCoeff_list[n];

            std::cout << "NumCoeff= " << NumCoeff<< std::endl;

            K_p_Coeff.resize(2*NumCoeff);
            K_p_ug.resize(2*NumCoeff);
            E_Coeff.resize(2*NumCoeff);
            DK_Coeff.resize(2*NumCoeff);
            DK_omega.resize(2*NumCoeff);
            N0.resize(2*NumCoeff);
            N0_tilde.resize(2*NumCoeff);
            N0_omega.resize(2*NumCoeff);
            htilde.resize(2*NumCoeff);
            Ttilde.resize(2*NumCoeff);
            BLtilde.resize(2*NumCoeff);
            BN.resize(2*NumCoeff);
            RBN.resize(2*NumCoeff);
            tempvector.resize(2*NumCoeff);
            Im_K.resize(2*NumCoeff);

            fill(K_p_Coeff.begin(), K_p_Coeff.end(), Vec2_int(0));
            fill(E_Coeff.begin(), E_Coeff.end(), Vec2_int(0));
            fill(DK_Coeff.begin(), DK_Coeff.end(), Vec2_int(0));
            fill(DK_omega.begin(), DK_omega.end(), Vec2_int(0));
            fill(N0.begin(), N0.end(), Vec2_int(0));
            fill(N0_tilde.begin(), N0_tilde.end(), Vec2_int(0));
            fill(N0_omega.begin(), N0_omega.end(), Vec2_int(0));
            fill(htilde.begin(), htilde.end(), 0);
            fill(Ttilde.begin(), Ttilde.end(), 0);
            fill(BLtilde.begin(), BLtilde.end(), 0);
            fill(BN.begin(), BN.end(), 0);
            fill(RBN.begin(), RBN.end(), 0);
            fill(tempvector.begin(), tempvector.end(), 0);
            fill(Im_K.begin(), Im_K.end(), 0);

            // ----------- Make Fourier coefficients into intervals ----------

            K_p_Coeff[0] = K_p_approx[0];
            K_p_Coeff[1] = K_p_approx[1];
            for (std::size_t k=1; k!=K_p_approx.size()/4; ++k) {
                K_p_Coeff[2*k] = K_p_approx[2*k];
                K_p_Coeff[2*k+1] = K_p_approx[2*k+1];
                K_p_Coeff[2*(NumCoeff-k)] = K_p_approx[2*(K_p_approx.size()/2-k)];
                K_p_Coeff[2*(NumCoeff-k)+1] = K_p_approx[2*(K_p_approx.size()/2-k)+1];
            }

            copy(K_p_Coeff.begin(),K_p_Coeff.end(),K_p_ug.begin());

            fourier_series_ug(K_p_ug);

            // Apply F
            for( std::size_t j=0; j!=NumCoeff; ++j) {
                INTERVAL temp = K_p_ug[2*j][0];
                Vec2_int z = K_p_ug[2*j];
                z[0] += INTERVAL(j)/INTERVAL(NumCoeff);
                //z[0] -= floor(z[0]);
                K_p_ug[2*j] = F_p(z,MP_int);
                K_p_ug[2*j][0] += temp - omega;
            }

            fourier_coeff_ug(K_p_ug);

            copy(K_p_Coeff.begin(),K_p_Coeff.end(),E_Coeff.begin());

            phase_translate_fourier(E_Coeff,omega);
            for (std::size_t k=0; k!=E_Coeff.size(); ++k)
                E_Coeff[k] -= K_p_ug[k];

            // Calculate P = [ DK N0 ]

            DK_Coeff[0][0] = 1;
            DK_Coeff[0][1] = 0;
            DK_Coeff[1][0] = 0;
            DK_Coeff[1][1] = 0;
            for (std::size_t k=1; k!=NumCoeff/2; ++k) {
                //std::cout << "k = " << k << std::endl;
                DK_Coeff[2*k] = -(twopi*k)*K_p_Coeff[2*k+1];
                DK_Coeff[2*k+1] = (twopi*k)*K_p_Coeff[2*k];

                DK_Coeff[2*(NumCoeff-k)] = (twopi*k)*K_p_Coeff[2*(NumCoeff-k)+1];
                DK_Coeff[2*(NumCoeff-k)+1] = -(twopi*k)*K_p_Coeff[2*(NumCoeff-k)];
            }
            DK_Coeff[NumCoeff] = -(twopi*NumCoeff)*K_p_Coeff[NumCoeff+1];
            DK_Coeff[NumCoeff+1] = (twopi*NumCoeff)*K_p_Coeff[NumCoeff];

            copy(DK_Coeff.begin(),DK_Coeff.end(),DK_omega.begin());
            phase_translate_fourier(DK_omega,omega);

            copy(DK_Coeff.begin(),DK_Coeff.end(),N0.begin());
            fourier_series_ug(N0);
            for (std::size_t j=0; j!=NumCoeff; ++j) {
                //std::cout << "j = " << j << std::endl;
                INTERVAL temp_int = N0[2*j][0]*N0[2*j][0] + N0[2*j][1]*N0[2*j][1];
                N0[2*j][0] /= temp_int;
                N0[2*j][1] /= temp_int;
                temp_int = N0[2*j][0];
                N0[2*j][0] = -N0[2*j][1];
                N0[2*j][1] = temp_int;
            }

            copy(DK_omega.begin(),DK_omega.end(),N0_omega.begin());
            fourier_series_ug(N0_omega);
            for (std::size_t j=0; j!=NumCoeff; ++j) {
                //std::cout << "j = " << j << std::endl;
                INTERVAL temp_int = N0_omega[2*j][0]*N0_omega[2*j][0] + N0_omega[2*j][1]*N0_omega[2*j][1];
                N0_omega[2*j][0] /= temp_int;
                N0_omega[2*j][1] /= temp_int;
                temp_int = N0_omega[2*j][0];
                N0_omega[2*j][0] = -N0_omega[2*j][1];
                N0_omega[2*j][1] = temp_int;
            }

            copy(N0.begin(),N0.end(),N0_tilde.begin());
            fourier_coeff_ug(N0_tilde);                 // OBS: THIS IS NOT RIGOROUS!!!
            //copy(N0_tilde.begin(),N0_tilde.end(),N0_omega.begin());
            //phase_translate_fourier(N0_omega,omega);

            fourier_series_ug(DK_omega);
            fourier_series_ug(N0_omega);

            copy(K_p_Coeff.begin(), K_p_Coeff.end(), K_p_ug.begin());
            fourier_series_ug(K_p_ug);

            for (std::size_t j=0; j!=NumCoeff; ++j) {
                INTERVAL epscos = MP_int.epsilon*cos(twopi*(INTERVAL(j)/NumCoeff+K_p_ug[2*j][0]));
                BLtilde[2*j] = DK_omega[2*j][0]-DK_omega[2*j][1];
                htilde[2*j] = (epscos*BLtilde[2*j]-DK_omega[2*j][1])*N0[2*j][0] + MP_int.lambda*(BLtilde[2*j]*N0[2*j][1]-1);
                BLtilde[2*j] = N0_omega[2*j][1] - N0_omega[2*j][0];
                Ttilde[2*j] = (epscos*BLtilde[2*j]-N0_omega[2*j][1])*N0[2*j][0] + MP_int.lambda*BLtilde[2*j]*N0[2*j][1];
            }

            for (std::size_t i=0; i!=2*NumCoeff; ++i)
                BN[i] = DK_Coeff[i][0] - DK_Coeff[i][1];
            phase_translate_fourier(BN,omega);

            // Calculate RBN = R_lambda BN
            
            copy(BN.begin(),BN.end(),RBN.begin());

            RBN[0] /= MP_int.lambda - 1;
            RBN[1] /= MP_int.lambda - 1;
            for (std::size_t k=1; k!=NumCoeff/2; ++k) {
                INTERVAL tempr, tempi, tempc, temp;
                tempi = sin(twopi*k*omega);
                tempr = cos(twopi*k*omega);
                tempc = MP_int.lambda*MP_int.lambda - 2*MP_int.lambda*tempr + 1;
                tempr = MP_int.lambda - tempr;
                temp = RBN[2*k];
                RBN[2*k] = tempr*temp - tempi*RBN[2*k+1];
                RBN[2*k+1] = tempi*temp + tempr*RBN[2*k+1];
                RBN[2*k] /= tempc;
                RBN[2*k+1] /= tempc;
                tempi = -tempi;
                temp = RBN[2*(NumCoeff-k)];
                RBN[2*(NumCoeff-k)] = tempr*temp - tempi*RBN[2*(NumCoeff-k)+1];
                RBN[2*(NumCoeff-k)+1] = tempi*temp + tempr*RBN[2*(NumCoeff-k)+1];
                RBN[2*(NumCoeff-k)] /= tempc;
                RBN[2*(NumCoeff-k)+1] /= tempc;
            }

            INTERVAL Phi_average(0);
            copy(RBN.begin(),RBN.end(),tempvector.begin());
            fourier_series_ug(tempvector);
            for( std::size_t j=0; j!=NumCoeff; ++j)
                Phi_average += BLtilde[2*j] - Ttilde[2*j]*tempvector[2*j];
            Phi_average = abs(Phi_average/NumCoeff);

            fourier_coeff_ug(htilde);
            fourier_coeff_ug(BLtilde);
            fourier_coeff_ug(Ttilde);

            for( std::size_t r=0; r!=rho_list.size() && !existence; ++r) {

                rho = rho_list[r];
                //rho_hat = 8*rho;
                rho_infinity = rho/4;
                delta = rho/4;
                std::cout << "rho = " << rho.right_l() << std::endl;
                std::cout << "rho_hat = " << rho_hat.right_l() << std::endl;
                std::cout << "rho_infinity = " << rho_infinity.right_l() << std::endl;
                std::cout << "delta = " << delta.right_l() << std::endl;

                INTERVAL C_NF = DFT_error_const(rho, rho_hat, NumCoeff);
                INTERVAL s_NF = DFT_error_coeff(0,rho_hat,NumCoeff);

                std::cout << "C_NF(rho_hat,rho) = " << C_NF.right_l() << std::endl;
                std::cout << "s_NF = " << s_NF.right_l() << std::endl;

                //INTERVAL DKrange = fourier_series_maxnorm_range_test(DK_Coeff, rho, 2*NumCoeff);
                //std::cout << "DKrange = " << DKrange << std::endl;
                INTERVAL DKrange = fourier_series_maxnorm_range(DK_Coeff, rho, 2*NumCoeff);
                //std::cout << "DKrange = " << DKrange << std::endl;

                INTERVAL B_L = DKrange.right_l();
                INTERVAL c_N = (INTERVAL(1)/DKrange).right_l();
                INTERVAL sigma_L = 2*B_L;
                INTERVAL B_P = 2*max(B_L,c_N);
                INTERVAL sigma_P = 2*B_P;

                DKrange = fourier_series_maxnorm_range(DK_Coeff, rho_hat, 2*NumCoeff);

                INTERVAL sigma_L_hat = DKrange.right_l();
                INTERVAL c_N_hat = (INTERVAL(1)/DKrange).right_l();

                INTERVAL normhtilde = fourier_max_norm(htilde,rho) + C_NF*((1+2*MP_int.epsilon+2*MP_int.lambda)*sigma_L_hat*c_N_hat + MP_int.lambda);

                INTERVAL B_H = 1/(1 - lambda - normhtilde);
                INTERVAL sigma_H = 2*B_H;

                INTERVAL e_T = C_NF*(1 + 2*epsilon + 2*lambda)*c_N_hat*c_N_hat;
                INTERVAL e_R = normhtilde/(1 - lambda)/(1 - lambda - normhtilde);
                INTERVAL e_BLtilde = 2*c_N*C_NF;

                if (e_R.right_l() < Real(0))
                    e_R = std::numeric_limits<double>::infinity();

                INTERVAL normBLtilde = fourier_max_norm(BLtilde,rho);
                INTERVAL normTtilde = fourier_max_norm(Ttilde,rho);
                INTERVAL normBN = fourier_max_norm(BN,rho);
                INTERVAL normRBN = fourier_max_norm(RBN,rho);

                std::cout << "normhtilde = " << normhtilde.right_l() << std::endl;
                std::cout << "e_T = " << e_T.right_l() << "\ne_R = " << e_R.right_l() << "\ne_BLtilde = " << e_BLtilde.right_l() << std::endl;
                std::cout << "normTtilde = " << normTtilde.right_l() << std::endl;
                std::cout << "normBN = " << normBN.right_l() << std::endl;
                std::cout << "normRBN = " << normRBN.right_l() << std::endl;

                INTERVAL r1star = normTtilde + e_T;
                INTERVAL r2star = Phi_average - s_NF*(normBLtilde + 2*c_N_hat*C_NF + r1star*normRBN);
                INTERVAL r3star = r1star*normBN*normhtilde/(1-lambda)/(1-lambda-normhtilde);

                INTERVAL B_D = 1/(r2star - r3star);
                INTERVAL sigma_D = 2*B_D;

                INTERVAL c_F1z = 1 + lambda + exp(twopi*rho), c_F1a = 1+INTERVAL(exp2(Real(-(int)prec))), c_F2 = twopi*exp(twopi*rho);

                INTERVAL C1 = 1 + sigma_H*sigma_P*c_N*c_F1z;

                INTERVAL C2 = 1 + sigma_P*sigma_D*c_F1a*C1;

                INTERVAL C3 = C2*((sigma_L + 1)*c_R*sigma_P*C1 + c_N*sigma_H*sigma_P*gamma*pow(delta,tau));

                INTERVAL C2_hat = sigma_L*C3 + c_N*sigma_H*sigma_P*C2*gamma*pow(delta,tau);

                INTERVAL C3_hat = sigma_D*sigma_P*C1;

                INTERVAL C23_hat = max(C2_hat,C3_hat*gamma*pow(delta,tau));

                INTERVAL C4_hat = 2*sigma_P*sigma_P*d*C2_hat;

                INTERVAL C4 = c_N*(sigma_P*c_F2*C23_hat*delta + c_F1z*C4_hat);

                INTERVAL C5_hat = 2*sigma_H*sigma_H*C4;

                INTERVAL C5 = C4_hat*c_F1a + sigma_P*c_F2*C23_hat;

                INTERVAL C6 = C5 + sigma_H*sigma_P*c_F1a*C4 + sigma_P*c_F1z*c_N*sigma_H*C5 + sigma_P*c_F1z*c_N*sigma_P*c_F1a*C5_hat;

                INTERVAL C6_hat = 2*sigma_D*sigma_D*C6;

                INTERVAL C7_hat = d*c_R*C3*gamma*pow(delta,tau-1) + c_F2*C23_hat/2;

                INTERVAL a1 = (rho-rho_infinity)/(rho-2*delta-rho_infinity), a3 = rho/delta;

                INTERVAL C8_hat = max( d*C2_hat/(sigma_L-B_L), C4_hat/(sigma_P-B_P) );
                C8_hat = max( C8_hat, C5_hat/(sigma_H-B_H) );
                C8_hat = max( C8_hat, C6_hat/(sigma_D-B_D) );
                C8_hat /= 1 - pow(a1,1-tau);
                C8_hat = max( C8_hat, C3_hat*pow(delta,tau+1)*gamma/max(mu,1-mu)/(1-pow(a1,2*tau)) );

                INTERVAL Cstar = max( C7_hat*pow(a1*a3,2*tau), C8_hat*pow(a3,tau+1)*gamma*pow(rho,tau-1) );

                INTERVAL Cdstar = max( C2_hat*pow(a3,tau)/(1 - pow(a1,-tau)), C3_hat*gamma*pow(rho,tau)/(1 - pow(a1,-2*tau)) );

                INTERVAL Ctilde = C23_hat*c_F2/2;

                INTERVAL Ctstar = pow(INTERVAL(4), tau)*Ctilde*Cdstar;

                for( std::size_t k=1; k!= NumCoeff_approx; ++k) {
                    INTERVAL temp = (exp(twopi*k*rho_hat)-exp(-twopi*k*rho_hat))/2;
                    Im_K[2*k] = -temp*K_p_Coeff[2*k+1][0];
                    Im_K[2*k+1] = temp*K_p_Coeff[2*k][0];
                    Im_K[2*(NumCoeff-k)] = temp*K_p_Coeff[2*(NumCoeff-k)+1][0];
                    Im_K[2*(NumCoeff-k)+1] = -temp*K_p_Coeff[2*(NumCoeff-k)][0];
                }
                INTERVAL rho_bar = fourier_series_real_range(Im_K,4*NumCoeff);
                rho_bar = max(-rho_bar.left_l(), rho_bar.right_l());

                INTERVAL FKbound = (1+MP_int.lambda)*sigma_L_hat + mu + epsilon*exp(twopi*(rho_hat+rho_bar))/twopi;

                INTERVAL normE = fourier_max_norm(E_Coeff,rho) + C_NF*FKbound;

                INTERVAL T1_const = Cstar*normE/(gamma*gamma*pow(rho,2*tau));

                INTERVAL T2_const = Cdstar*normE/(gamma*pow(rho,tau));

                INTERVAL T3_const1 = Ctstar*normE/(gamma*gamma*pow(rho_infinity,tau)*pow(rho,tau));

                INTERVAL T3_const2 = Cdstar*(gamma*pow(rho_infinity,tau)/Ctstar - normE/(gamma*pow(rho,tau)));

                std::cout << "sigma_L = " << sigma_L.right_l() << "\nc_N = " << c_N.right_l() << "\nsigma_P = " << sigma_P.right_l() << std::endl;
                std::cout << "sigma_H = " << sigma_H.right_l() << std::endl;
                //std::cout << "E_D = " << E_D.right_l() << std::endl;
                std::cout << "sigma_D = " << sigma_D.right_l() << std::endl;
                std::cout << "rho_bar = " << rho_bar.right_l() << std::endl;
                std::cout << "FKbound = " << FKbound.right_l() << std::endl;
                std::cout << "normE = " << normE.right_l() << std::endl;

                std::cout << "T1_const = " << T1_const.right_l() << "\n";

                if (T1_const.right_l() < my_mpfr(1)) {
                    existence = true;
                    std::cout << "\nThere exists an invariant circle\n";
                    std::cout << "Distance to approximation is less than " << T2_const.left_l() << std::endl;
                    if (T3_const1.right_l() < my_mpfr(1))
                        std::cout << "Local uniquness within radius " << T3_const2.left_l() << "\n" << std::endl;
                } else
                    std::cout << "\nCould not prove existance of invariant circle\n" << std::endl;

            }
        }




    }



    return 0;
}


