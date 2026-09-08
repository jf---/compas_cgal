// Does the filter's value scale with EXPRESSION DEPTH / COEFFICIENT GROWTH?
// Chain rational constructions (as a fitted centre or projection chain does),
// then make ONE decision on the result. Sweep the chain depth.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/CORE/BigRat.h>
#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;
using LazyFT = Epeck::FT;

template <typename NT>
NT chain(NT x, const NT& a, const NT& b, int depth) {
    for (int i = 0; i < depth; ++i) x = (x * a + b) / (x + a);   // rational, coefficients grow
    return x;
}

template <typename NT, typename Conv>
double run(const std::vector<double>& seeds, int depth, Conv conv, int& sink) {
    const NT a = conv(3.7182818), b = conv(-1.4142136), zero = conv(0.0);
    std::vector<NT> s;
    s.reserve(seeds.size());
    for (double d : seeds) s.push_back(conv(d));
    const auto t0 = std::chrono::steady_clock::now();
    for (const NT& x0 : s) {
        const NT v = chain<NT>(x0, a, b, depth);
        if (v > zero) ++sink; else --sink;            // ONE generic decision
    }
    const auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(t1 - t0).count();
}

int main() {
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> dist(-100.0, 100.0);
    std::vector<double> seeds(200);
    for (double& d : seeds) d = dist(rng);
    std::printf("%6s %14s %14s %10s\n", "depth", "lazy us/dec", "BigRat us/dec", "ratio");
    for (int depth : {1, 2, 4, 8, 12, 16}) {
        int sink = 0;
        const double tl = run<LazyFT>(seeds, depth, [](double d){ return LazyFT(d); }, sink);
        const double tc = run<CORE::BigRat>(seeds, depth, [](double d){ return CORE::BigRat(d); }, sink);
        std::printf("%6d %14.2f %14.2f %9.2fx\n", depth,
                    1e6 * tl / seeds.size(), 1e6 * tc / seeds.size(), tc / tl);
    }
    return 0;
}
