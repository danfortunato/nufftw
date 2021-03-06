# NUFFTW
The **N**on**U**niform **F**astest **F**ourier **T**ransform in the **W**est

NUFFTW is a C++ library for computing the nonuniform discrete Fourier transform (NUDFT) of types I, II, & III in O(N log N) time, where N is the size of the transform. It is based on the paper:

* [1] D. Ruiz-Antolin and A. Townsend, "A nonuniform fast Fourier transform based on low rank approximation," submitted, 2017.

This paper relates the NUFFT to sums of FFTs by low rank approximation.
