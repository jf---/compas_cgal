// Question #2: which point-location strategy fits the audit's mutate/query
// interleave? Pattern = (insert 1 curve; do Q locate() queries) x N.
// For each strategy measure (a) total INSERT time with the PL attached -- reveals
// whether it rebuilds or updates incrementally on each change -- and (b) total
// QUERY time. Segment traits: PL rebuild/incremental behavior is traits-agnostic.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <chrono>
#include <iostream>
#include <random>

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits = CGAL::Arr_segment_traits_2<K>;
using Arr = CGAL::Arrangement_2<Traits>;
using Point = Traits::Point_2;

template <class PL>
void bench(const char* name, int N, int Q) {
    Arr arr;
    PL pl(arr);  // dynamic strategies register as observers here
    Traits traits;
    auto mk = traits.construct_x_monotone_curve_2_object();
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> U(0.0, 1000.0);

    double insert_ms = 0, query_ms = 0;
    for (int i = 0; i < N; ++i) {
        double x0 = U(rng), y0 = U(rng), x1 = U(rng), y1 = U(rng);
        if (x0 == x1 && y0 == y1) x1 += 1.0;
        auto t0 = std::chrono::steady_clock::now();
        CGAL::insert(arr, mk(Point(x0, y0), Point(x1, y1)));
        auto t1 = std::chrono::steady_clock::now();
        insert_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();

        auto q0 = std::chrono::steady_clock::now();
        for (int j = 0; j < Q; ++j) {
            auto r = pl.locate(Point(U(rng), U(rng)));
            (void)r;
        }
        auto q1 = std::chrono::steady_clock::now();
        query_ms += std::chrono::duration<double, std::milli>(q1 - q0).count();
    }
    std::cout << name << ": edges=" << arr.number_of_edges()
              << "  insert_total=" << insert_ms << "ms"
              << "  query_total=" << query_ms << "ms\n";
}

int main() {
    int N = 250, Q = 20;
    std::cout << "pattern: (insert 1 seg; " << Q << " locate) x " << N << "\n";
    bench<CGAL::Arr_naive_point_location<Arr>>("naive        ", N, Q);
    bench<CGAL::Arr_walk_along_line_point_location<Arr>>("walk         ", N, Q);
    // landmarks: CGAL assertion crash (new_face != face) under heavy dynamic
    // insertion -- unsuitable for a constantly-mutating arrangement. Skipped.
    bench<CGAL::Arr_trapezoid_ric_point_location<Arr>>("trapezoid_ric", N, Q);
    return 0;
}
