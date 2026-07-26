#include "segment_site_mat.h"

#include <cstdlib>

int main()
{
    const SegmentSiteMatCompileEvidence2 evidence =
        segment_site_mat_compile_spike();
    return evidence.delaunay_valid
            && evidence.assigned_dual_primitives > 0
            && evidence.exact_clearance_roots == 2
        ? EXIT_SUCCESS
        : EXIT_FAILURE;
}
