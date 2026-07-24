// Minimal compile-test: instantiate an Arrangement_2, insert one segment.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Arr_segment_traits_2<K> Traits;
typedef CGAL::Arrangement_2<Traits> Arr;

int main() {
    Arr arr;
    Traits::Point_2 p(0, 0), q(1, 1);
    Traits traits;
    auto mk = traits.construct_x_monotone_curve_2_object();
    CGAL::insert(arr, mk(p, q));
    std::cout << "OK faces=" << arr.number_of_faces()
              << " edges=" << arr.number_of_edges()
              << " vertices=" << arr.number_of_vertices() << "\n";
    return 0;
}
