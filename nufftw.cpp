#include <cmath>
//#include <tuple>
//#include <armadillo>
#include <gsl/gsl_sf_bessel.h>  // Bessel functions
#include <gsl/gsl_sf_lambert.h> // Lambert W functions
#include "nufftw.h"

namespace nufftw
{
    typedef arma::Mat<double>   Mat;
    typedef arma::Mat<complex> CMat;
    typedef arma::Col<double>   Vec;
    typedef arma::Col<complex> CVec;

    /** Evaluate Chebyshev polynomials of degree 0,...,n at x. */
    Mat chebT(int n, const Vec& x)
    {
        int N = x.n_rows;
        Mat T = arma::zeros(N,n+1);
        T.col(0).ones();
        if (n>0) {
            T.col(1) = x;
            for (int k=1; k<n; ++k) {
                T.col(k+1) = 2*x%T.col(k)-T.col(k-1);
            }
        }
        return T;
    }

    /** The bivarate Chebyshev coefficients for the function
     *  f(x,y) = exp(-i*x.*y) on the domain [-gam, gam]x[0,2*pi], as given by
     *  Lemma A.2 of Townsend's DPhil thesis. */
    void bessel_coeffs(int K, double gamma, complex* cfs)
    {
        const complex I(0,1);
        //CMat cfs(K,K,arma::fill::zeros);
        double arg = -0.5*gamma*M_PI;
        for (int p=0; p<K; ++p) {
            for (int q=p%2; q<K; q+=2) {
                cfs(p,q) = 4. * std::pow(I,q) * gsl_sf_bessel_Jn((p+q)*0.5, arg) * gsl_sf_bessel_Jn((q-p)*0.5, arg);
            }
        }
        cfs.row(0)*=0.5;
        cfs.col(0)*=0.5;
        return cfs;
    }

    // arg = -gam*pi/2;
    // [pp, qq] = meshgrid(0:K-1);
    // cfs = 4*(1i).^qq.*besselj((pp+qq)/2,arg).*besselj((qq-pp)/2, arg);
    // cfs(2:2:end,1:2:end) = 0;
    // cfs(1:2:end,2:2:end) = 0;
    // cfs(1,:) = cfs(1,:)/2;
    // cfs(:,1) = cfs(:,1)/2;

    // std::tuple<Vec, Vec, double, int> parameters(const Vec& x, double tol)
    // {
    //     int N = x.n_rows;
    //     Vec s = arma::round(N*x);
    //     Vec t = s - arma::floor(s/N)*N;
    //     double gamma = arma::norm(N*x - s, "inf");
    //     double xi = std::log(std::log(10/tol)/(7*gamma));
    //     double lw = gsl_sf_lambert_W0(xi);
    //     int K = std::ceil(5*gamma*std::exp(lw));
    //     // double log_xi = std::log(xi);
    //     // double inv_xi = 1.0/xi;
    //     // double lw = xi - log_xi*(1.0 + inv_xi*(1.0 + inv_xi*(0.5*log_xi - 1.0)));
    //     return {s, t, gamma, K};
    // }

    Vec closestGridpoint(const Vec& x)
    {
        int N = x.n_rows;
        return arma::round(N*x);
    }

    Vec closestFFTSample(const Vec& s)
    {
        int N = s.n_rows;
        return s - arma::floor(s/N)*N;
    }

    double perturbationParameter(const Vec& x, const Vec& s)
    {
        int N = x.n_rows;
        return arma::norm(N*x-s, "inf");
    }

    int rankParameter(double gamma, double tol)
    {
        double xi = std::log(std::log(10./tol)/(7.*gamma));
        double lw = gsl_sf_lambert_W0(xi);
        return std::ceil(5.*gamma*std::exp(lw));
    }

    /** Construct a low rank approximation to A_{jk} = exp(-2*pi*1i*(x_j-s_j/N)*k),
     *  0<=j,k<=N-1, where |x_j-j/N|<=gamma<=1/2. */
    // std::tuple<Mat<complex>, Mat<complex>> constructUV(const Vec& x, const Vec& omega, int K, double tol)
    // {
    //     int N = omega.n_rows;
    //     auto [s, t, gamma, K] = parameters(x, tol);
    //     Vec er = N*x - s;
    //     CVec scl = arma::exp(-I*M_PI*er);
    //     CMat u = (chebT(K-1, er/gamma) * besselCoeffs(K, gamma)).each_col() % scl;
    //     CMat v = chebT(K-1, 2*omega/N-1);
    // }

    CMat constructU(const Vec& x, const Vec& s, int gamma,  int K)
    {
        const complex I(0,1);
        int N = x.n_rows;
        Vec er = N*x - s;
        CVec scl = arma::exp(-I*M_PI*er);
        CMat u = (chebT(K-1, er/gamma) * besselCoeffs(K, gamma)).each_col() % scl;
        return u;
    }

    CMat constructV(int N, int K)
    {
        CMat v = chebT(K-1, 2*arma::linspace(0,N-1)/N-1);
        return v;
    }

    CMat constructV(const Vec& omega, int K)
    {
        int N = omega.n_rows;
        CMat v = chebT(K-1, 2*omega/N-1);
        return v;
    }

    plan_nufft2_1d()
    nufftw::execute();
    nufftw::destroy_plan()

    nufftw::plan_nudft2_1d

    nufftw_plan nufft2(int N, const complex* c, const double* x, complex* f, double tol)
    {
        // Compute parameters
        Vec s = closestGridpoint(x);
        Vec t = closestFFTSample(s);
        double gamma = perturbationParameter(x, s);
        int K = rankParameter(gamma, tol);

        CMat u = constructU(x, s, gamma, K);
        CMat v = constructV(N, K);

        // nufftw_plan p = plan_nufft2(N, x);
        nufftw_plan p = [&]() {

        };
        p(c,f);
        return p;
    }
}
