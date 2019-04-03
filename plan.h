#ifndef NUFFTW_PLAN_H
#define NUFFTW_PLAN_H

#include <complex>
#include <cmath>
#include <limits>               // st::numeric_limits
#include <algorithm>            // std::max
#include <gsl/gsl_sf_bessel.h>  // Bessel functions
#include <gsl/gsl_sf_lambert.h> // Lambert W functions
#include <fftw3.h>
#include "timer.h"

namespace nufftw
{
    using namespace std::complex_literals;
    typedef std::complex<double> complex;

    struct options
    {
        double tol = std::numeric_limits<double>::epsilon();
        unsigned int fftw_flags = FFTW_ESTIMATE;
    };

    namespace details
    {
        void cheb_vals(int n, int p, const double* x, double* vals);
        void bessel_coeffs(int r, double gam, complex* cfs);
    }

    class plan_base
    {
    	public:
    		plan_base(const int* ns_, complex* in_, complex* out_, const options& opts_ = options()) :
                in(in_), out(out_), opts(opts_)
            {
                n=1;
                for (int i=0; i<N; ++i) {
                    ns[i]=ns_[i];
                    n*=ns[i];
                }
            }

            virtual ~plan_base()
            {
                fftw_destroy_plan(fp);
            }

    		virtual void execute() = 0;

        protected:
            int m;               // Number of nonuniform points
            int ns[N];           // Size of the problem in each dimension
            int n;               // Total size of the problem
            complex* in;         // Input coefficients
            complex* out;        // Output values
            fftw_plan fp;        // FFTW plans
            options opts;        // Options
            double gam;          // Perturbation parameter
            int r;               // Rank parameter
    };

    class plan1 : public plan_base<N>
    {
        public:
            plan1(const int* n_, complex* in_, complex* out_, const double* omega_, const options& opts_ = options()) :
                plan_base<N>(n_, in_, out_, opts_), omega(omega_)
            {}
            virtual ~plan1();
            virtual void execute();
        protected:
            const double* omega;
            plan1<1>* plans;
    };

    class plan1_2d : public plan_base
    {
        
    };

    class plan1_1d : public plan_base
    {
        public:
            plan1(int n_, complex* in_, complex* out_, const double* omega_, const options& opts_ = options()) :
                plan_base<1>(&n_, in_, out_, opts_), omega(omega_), t(new int[n]),
                er(new double[n]), temp(new double[n])
            {
                timer::tic();
                gam = 0;
                for (int i=0; i<n; ++i) {
                    int s = std::round(omega[i]);
                    t[i] = s%n;
                    er[i] = omega[i]-s;
                    gam = std::max(gam, std::abs(er[i]));
                }
                timer::toc("Compute gamma");

                timer::tic();
                r = 1;
                if (gam>opts.tol) {
                    // Asymptotic approximation to Lambert-W:
                    //   double xi = std::log(std::log(10./opts.tol)/(7.*gam);
                    //   double log_xi=std::log(xi), inv_xi=1./xi;
                    //   double lw = xi-log_xi*(1+inv_xi*(1+inv_xi*(0.5*log_xi-1)));
                    double lw = gsl_sf_lambert_W0(std::log(10./opts.tol)/(7.*gam));
                    r = std::ceil(5.*gam*std::exp(lw));
                }
                timer::toc("Compute r");

                timer::tic();
                // Plan the FFTs
                fft_data = reinterpret_cast<complex*>(fftw_alloc_complex(n*r));
                // int stride = r;
                int stride = 1;
                // int dist = 1;
                int dist = n;
                fp = fftw_plan_many_dft(1, &n, r,
                    reinterpret_cast<fftw_complex*>(fft_data), NULL, stride, dist,
                    reinterpret_cast<fftw_complex*>(fft_data), NULL, stride, dist,
                    FFTW_BACKWARD, opts.fftw_flags
                );
                timer::toc("Plan the FFTs");

                // Construct a low rank approximation, using Chebyshev expansions
                // for A_K = exp(-2*pi*1im*(x[j]-j/N)*k):
                
                // Construct a low rank approximation to
                //     A_{jk} = exp(-2*pi*1i*(x_j-s_j/N)*k), 0<=j,k<=N-1,
                // where |x_j-j/N|<= gam <=1/2.  See [1].
                u = new complex[n*r];
                v = new double[n*r];

                // Compute u
                timer::tic();
                double* cheb = new double[n*r];
                complex* bess = new complex[r*r];
                timer::tic();
                double fac = 1./gam;
                for (int i=0; i<n; ++i) temp[i]=er[i]*fac;
                details::cheb_vals(n, r, temp, cheb);
                timer::toc("Chebyshev");
                timer::tic();
                details::bessel_coeffs(r, gam, bess);
                timer::toc("Bessel");
                timer::tic();
                for (int i=0; i<n; ++i) {
                    complex scl = std::exp(-1i*M_PI*er[i]);
                    for (int j=0; j<r; ++j) {
                        u[r*i+j]=0;
                        for (int k=0; k<r; ++k) {
                            u[r*i+j]+=cheb[r*i+k]*bess[r*j+k];
                        }
                        u[r*i+j]*=scl;
                    }
                }
                timer::toc("Matrix multiply");
                delete [] bess;
                delete [] cheb;
                timer::toc("Compute u");

                // Compute v
                timer::tic();
                fac = 2./n;
                temp[0]=-1;
                for (int i=1; i<n; ++i) temp[i]=temp[i-1]+fac;
                details::cheb_vals(n, r, temp, v);
                timer::toc("Compute v");
            }

            virtual ~plan1()
            {
                delete [] t;
                delete [] er;
                delete [] u;
                delete [] v;
                delete [] temp;
                fftw_free(reinterpret_cast<fftw_complex*>(fft_data));
            }

            virtual void execute()
            {
                for (int i=0; i<n*r; ++i) fft_data[i]=0;
                for (int i=0; i<n; ++i) {
                    for (int j=0; j<r; ++j) {
                        // fft_data[r*t[i]+j] += std::conj(in[i]*u[r*i+j]);
                        fft_data[n*j+t[i]] += std::conj(in[i]*u[r*i+j]);
                    }
                }

                // Do the FFTs
                timer::tic();
                fftw_execute(fp); // normalize?
                timer::toc("Call FFTW");

                for (int i=0; i<n; ++i) {
                    out[i]=0;
                    for (int j=0; j<r; ++j) {
                        // out[i] += v[r*i+j]*std::conj(fft_data[r*i+j]);
                        out[i] += v[r*i+j]*std::conj(fft_data[n*j+i]);
                    }
                }
            }

        protected:
            const double* omega; // FFT samples
            int* t;              // Closest equispaced FFT samples to omega
            double* er;          // Perturbation error
            complex* u;          // Low rank approximation
            double* v;           // Low rank approximation
            double* temp;        // Temporary storage
            complex* fft_data;   // Data for FFTW
    };

    // template<int N>
    // class plan2 : public plan_base<N>
    // {
    //     public:
    //         plan2();
    //         virtual ~plan2();
    // };

    // template<>
    // class plan2<1> : public plan_base<1>
    // {
    //     public:
    //         plan2();
    //         virtual ~plan2();
    // };

    // template<int N>
    // class plan3 : public plan_base<N>
    // {
    //     public:
    //         plan3();
    //         virtual ~plan3();
    // };

    // template<>
    // class plan3<1> : public plan_base<1>
    // {
    //     public:
    //         plan3();
    //         virtual ~plan3();
    // };

    typedef plan_base* plan;

    // Methods for planning
    template<int N>
    void destroy_plan(plan<N> p)
    {
        if (p) delete p;
    }

    template<int N>
    void execute(const plan<N> p)
    {
        p->execute();
    }

    namespace details
    {
        /** Evaluate Chebyshev polynomials of degree 0,...,p-1 at x.
         *  (row-major order) */
        void cheb_vals(int n, int p, const double* x, double* vals)
        {
            for (int i=0; i<n; ++i) {
                vals[p*i]=1;
                vals[p*i+1]=x[i];
            }
            for (int i=0; i<n; ++i) {
                for (int j=1; j<p-1; ++j) {
                    vals[p*i+(j+1)] = 2.*x[i]*vals[p*i+j] - vals[p*i+(j-1)];
                }
            }
        }

        /** The bivarate Chebyshev coefficients for the function
         *  f(x,y) = exp(-i*x.*y) on the domain [-gam, gam]x[0,2*pi], as given
         *  by Lemma A.2 of Townsend's DPhil thesis. (column-major order) */
        void bessel_coeffs(int r, double gam, complex* cfs)
        {
            double arg = -0.5*gam*M_PI;
            for (int p=0; p<r; ++p) {
                for (int q=p%2; q<r; q+=2) {
                    cfs[r*q+p] = 4.*std::pow(1i,q)*gsl_sf_bessel_Jn((p+q)*0.5,arg)*gsl_sf_bessel_Jn((q-p)*0.5,arg);
                }
            }
            for (int p=0; p<r; ++p) cfs[p]*=0.5;
            for (int q=0; q<r; ++q) cfs[r*q]*=0.5;
        }
    }
}

#endif
