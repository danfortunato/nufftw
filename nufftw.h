#ifndef NUFFTW_H
#define NUFFTW_H

#include <complex>
#include <functional>
#include <limits>
#include <fftw3.h>

#include "plan.h"

namespace nufftw
{
    // typedef std::complex<double> complex;
    // typedef std::function<void(const double*, double*)> nufftw_plan;
    // plan nufft1(int n, const complex* in, complex* out, const double* omega, double tol = std::numeric_limits<double>::epsilon());
    // plan nufft2(int n, const complex* in, complex* out, const double* x, double tol = std::numeric_limits<double>::epsilon());
    // plan nufft3(int n, const complex* in, complex* out, const double* x, const double* omega, double tol = std::numeric_limits<double>::epsilon());

    /* NUFFT-I */

    // Direct
    //template<int N>
    //plan<N> nufft1(const int* n,              const complex* in, complex* out, const double* omega, const options& opts = options());
    // plan<1> nufft1_1d(int n0, complex* in, complex* out, const double* omega, const options& opts = options())
    // {
    //     plan<1> p = plan_nufft1_1d(n0, in, out, omega, opts);
    //     execute(p);
    //     return p;
    // }
    //plan<2> nufft1_2d(int n0, int n1,         const complex* in, complex* out, const double* omega, const options& opts = options());
    //plan<3> nufft1_3d(int n0, int n1, int n2, const complex* in, complex* out, const double* omega, const options& opts = options());
    // Planners
    //template<int N>
    //plan<N> plan_nufft1(const int* n,              const complex* in, complex* out, const double* omega, const options& opts = options());
    plan<1> plan_nufft1_1d(int n0, complex* in, complex* out, const double* omega, const options& opts = options())
    {
        return new plan1<1>(n0, in, out, omega, opts);
    }
    void nufft1_1d(int n0, complex* in, complex* out, const double* omega, const options& opts = options())
    {
        plan<1> p = plan_nufft1_1d(n0, in, out, omega, opts);
        execute(p);
        destroy_plan(p);
    }
    plan<1> plan_many_nufft1_1d(int n0, int howmany,
                             complex *in,  int inembed, int istride, int idist,
                             complex *out, int onembed, int ostride, int odist,
                             const double* omega, const options& opts = options())
    {
        return new plan1<1>(n0, howmany,
            in,  inembed, istride, idist,
            out, onembed, ostride, odist,
            omega, opts);
    }
    //plan<2> plan_nufft1_2d(int n0, int n1,         const complex* in, complex* out, const double* omega, const options& opts = options());
    //plan<3> plan_nufft1_3d(int n0, int n1, int n2, const complex* in, complex* out, const double* omega, const options& opts = options());

    // /* NUFFT-II */

    // // Direct
    // template<int N>
    // plan<N> nufft2(const int* n,              const complex* in, complex* out, const double* x, const options& opts = options());
    // plan<1> nufft2_1d(int n0,                 const complex* in, complex* out, const double* x, const options& opts = options());
    // plan<2> nufft2_2d(int n0, int n1,         const complex* in, complex* out, const double* x, const options& opts = options());
    // plan<3> nufft2_3d(int n0, int n1, int n2, const complex* in, complex* out, const double* x, const options& opts = options());
    // // Planners
    // template<int N>
    // plan<N> plan_nufft2(const int* n,              const complex* in, complex* out, const double* x, const options& opts = options());
    // plan<1> plan_nufft2_1d(int n0,                 const complex* in, complex* out, const double* x, const options& opts = options());
    // plan<2> plan_nufft2_2d(int n0, int n1,         const complex* in, complex* out, const double* x, const options& opts = options());
    // plan<3> plan_nufft2_3d(int n0, int n1, int n2, const complex* in, complex* out, const double* x, const options& opts = options());

    // /* NUFFT-III */

    // // Direct
    // template<int N>
    // plan<N> nufft3(const int* n,              const complex* in, complex* out, const double* x, const double* omega, const options& opts = options());
    // plan<1> nufft3_1d(int n0,                 const complex* in, complex* out, const double* x, const double* omega, const options& opts = options());
    // plan<2> nufft3_2d(int n0, int n1,         const complex* in, complex* out, const double* x, const double* omega, const options& opts = options());
    // plan<3> nufft3_3d(int n0, int n1, int n2, const complex* in, complex* out, const double* x, const double* omega, const options& opts = options());
    // // Planners
    // template<int N>
    // plan<N> plan_nufft3(const int* n,              const complex* in, complex* out, const double* x, const double* omega, const options& opts = options());
    // plan<1> plan_nufft3_1d(int n0,                 const complex* in, complex* out, const double* x, const double* omega, const options& opts = options());
    // plan<2> plan_nufft3_2d(int n0, int n1,         const complex* in, complex* out, const double* x, const double* omega, const options& opts = options());
    // plan<3> plan_nufft3_3d(int n0, int n1, int n2, const complex* in, complex* out, const double* x, const double* omega, const options& opts = options());
}

#endif
