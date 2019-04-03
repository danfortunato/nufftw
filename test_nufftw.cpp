#include <complex>
#include <random>
#include <iostream>
#include "nufftw.h"
#include "timer.h"
using namespace std::complex_literals;

typedef std::complex<double> complex;

// Exact solutions
void nudft1(int n, const complex* in, complex* out, const double* omega);
// void nudft2(int N, const complex* c, const real* x,                    complex* f);
// void nudft3(int N, const complex* c, const real* x, const real* omega, complex* f);

int main()
{
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<> dist(0,1);

    nufftw::options opts;
    opts.tol = 1e-5;
    int n = 10000000;
    complex* exact = new complex[n];
    complex* fast  = new complex[n];
    complex* c     = new complex[n];
    double* omega  = new double[n];
    
    for (int i=0; i<n; ++i) {
        omega[i] = i + 1.2*dist(rng)/n;
        c[i] = complex(dist(rng),dist(rng));
    }

    // nudft1(n, c, exact, omega);

    timer::tic();
    timer::tic();
    auto p = nufftw::plan_nufft1_1d(n, c, fast, omega, opts);
    timer::toc("Plan NUFFT");

    timer::tic();
    nufftw::execute(p);
    timer::toc("Execute NUFFT");
    timer::toc();

    nufftw::destroy_plan(p);

    // double norm=0;
    // for (int i=0; i<n; ++i) {
    //     double err = std::abs(exact[i]-fast[i]);
    //     if (err>norm) norm=err;
    // }
    // std::cout << norm << std::endl;

    delete [] c;
    delete [] omega;
    delete [] fast;
    delete [] exact;

    // double tol = 1e-15; // The desired tolerance
    // arma_rng::set_seed_random();
    // std::vector<bool> pass;

    // // Test NUFFT-II
    // int N=1;
    // for (int i=0; i<5; ++i, N*=10) {
    //     Vec x = linspace<Vec>(0,1,N+1);
    //     x.shed_row(N);
    //     x+=1.2/N;
    //     CVec c(N), exact(N), fast(N);
    //     c.randu();
    //     nudft2(N, c.memptr(), x.memptr(), exact.memptr());
    //     nufft2(N, c.memptr(), x.memptr(), fast.memptr());
    //     pass.push_back(norm(exact-fast, "inf") < 300*N*tol*arma::norm(c,1));
    // }

    // return 0;
}

void nudft1(int n, const complex* in, complex* out, const double* omega)
{
    double fac=1./n;
    for (int j=0; j<n; ++j) {
        out[j]=0;
        double jn=j*fac;
        for (int k=0; k<n; ++k) {
            out[j] += std::exp(-2*M_PI*1i*jn*omega[k])*in[k];
        }
    }
}

// void nudft2(int N, const complex* c, const real* x, complex* f)
// {
//     for (int j=0; j<N; ++j) {
//         for (int k=0; k<N; ++k) {
//             f[j] = std::exp(-2*M_PI*1i*x[j]*k)*c[k];
//         }
//     }
// }

// void nudft3(int N, const complex* c, const real* x, const real* omega, complex* f)
// {
//     for (int j=0; j<N; ++j) {
//         for (int k=0; k<N; ++k) {
//             f[j] = std::exp(-2*M_PI*1i*x[j]*omega[k])*c[k];
//         }
//     }
// }
