#include <vector>
#include <chrono>
#include <stdio.h>
#include <fftw3.h>

int main()
{
    int ntrials = 10;
    std::vector<int> ns = {4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152};
    std::vector<std::chrono::duration<double>> ts;

    for (int n : ns) {
        auto start = std::chrono::high_resolution_clock::now();
        for (int i=0; i<ntrials; ++i) {
            fftw_complex *in, *out;
            fftw_plan p;
            in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
            out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
            p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(p);
            fftw_destroy_plan(p);
            fftw_free(in);
            fftw_free(out);
        }
        auto end = std::chrono::high_resolution_clock::now();
        ts.push_back(end-start);
        printf("# n=%d: %gs\n", n, ts.back().count()/ntrials);
    }
}
