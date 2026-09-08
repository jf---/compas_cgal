// Measures the ONE claim the remediation rests on:
// the same decision, on an unfiltered exact rational vs a lazy/filtered one.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/CORE/BigRat.h>
#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;
using LazyFT = Epeck::FT;                       // Lazy_exact_nt<cpp_rational>

// A generic orientation-shaped decision: sign of a 2x2 determinant of differences,
// then a squared-distance comparison. Both are ordinary deciding work.
template <typename NT>
int decide(const NT& ax, const NT& ay, const NT& bx, const NT& by,
           const NT& cx, const NT& cy, const NT& r2) {
    const NT d = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
    const NT s = (cx - ax) * (cx - ax) + (cy - ay) * (cy - ay);
    int v = 0;
    if (d > NT(0)) v += 1; else if (d < NT(0)) v -= 1;
    if (s > r2) v += 2; else if (s < r2) v -= 2;
    return v;
}

template <typename NT, typename Conv>
double run(const std::vector<double>& coords, int reps, Conv conv, int& sink) {
    std::vector<NT> c;
    c.reserve(coords.size());
    for (double d : coords) c.push_back(conv(d));
    const NT r2 = conv(1234.5678);
    const auto t0 = std::chrono::steady_clock::now();
    for (int k = 0; k < reps; ++k)
        for (size_t i = 0; i + 6 < c.size(); i += 6)
            sink += decide<NT>(c[i], c[i+1], c[i+2], c[i+3], c[i+4], c[i+5], r2);
    const auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(t1 - t0).count();
}

int main() {
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(-100.0, 100.0);
    std::vector<double> coords(6 * 2000);
    for (double& d : coords) d = dist(rng);           // GENERIC doubles, not integers
    const int reps = 20;
    const long decisions = (long)reps * (coords.size() / 6) * 2;

    int sink = 0;
    const double t_lazy = run<LazyFT>(coords, reps, [](double d){ return LazyFT(d); }, sink);
    const double t_core = run<CORE::BigRat>(coords, reps, [](double d){ return CORE::BigRat(d); }, sink);

    std::printf("decisions            %ld\n", decisions);
    std::printf("Lazy_exact_nt<cpp_rational>  %8.4f s   %9.1f ns/decision\n",
                t_lazy, 1e9 * t_lazy / decisions);
    std::printf("bare CORE::BigRat            %8.4f s   %9.1f ns/decision\n",
                t_core, 1e9 * t_core / decisions);
    std::printf("ratio (BigRat / lazy)        %8.2fx\n", t_core / t_lazy);
    std::printf("sink %d\n", sink);
    return 0;
}
