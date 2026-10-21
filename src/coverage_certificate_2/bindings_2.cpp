// nanobind module `_coverage_certificate_2`: the CORE-free coverage certificate.
//
// Exact rationals cross this boundary as (signed decimal numerator, positive
// decimal denominator) text pairs, `Fraction(int(n), int(d))` in Python. No
// double leaves this boundary except the reporting `elapsed_s`.

#include "coverage_certificate_2/certificate_2.h"
#include "coverage_certificate_2/circle_sweep_2.h"
#include "coverage_certificate_2/errors_2.h"
#include "coverage_certificate_2/exact_certificate_2.h"
#include "coverage_certificate_2/line_sweep_2.h"
#include "coverage_certificate_2/tower_2.h"
#include "exact/errors.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace nb = nanobind;
using namespace nb::literals;
namespace cc = compas_cgal::coverage_certificate;

namespace {
using XY = std::array<double, 2>;
using RationalPair = std::pair<std::string, std::string>;
using PointText = std::pair<RationalPair, RationalPair>;
using EdgeRecord = std::tuple<std::uint32_t, std::uint32_t, int, std::uint32_t, bool, std::vector<std::uint32_t>,
                              std::vector<std::uint32_t>, RationalPair, PointText, PointText, bool>;
using WitnessRecord = std::tuple<int, std::uint32_t, std::uint32_t, std::uint32_t, RationalPair, PointText>;
/// (is_line, sense, values): a declared primitive; values are
/// line: (sx, sy, ex, ey, ux, uy, length), arc: (sx, sy, ex, ey, cx, cy, radius).
using PrimitiveRow = std::tuple<bool, int, std::vector<RationalPair>>;
/// (is_line, sense, support, source, target): a dOmega* piece; support is
/// line: (nx, ny, level), circle: (cx, cy, R, rho); source/target are one-root
/// points (x0, x1, y0, y1, radicand).
using PieceRow = std::tuple<bool, int, std::vector<RationalPair>, std::vector<RationalPair>, std::vector<RationalPair>>;
/// (loop, segment0, tau0, segment1, tau1): a connector forward along one loop.
using ConnectorRow = std::tuple<std::uint32_t, std::uint32_t, RationalPair, std::uint32_t, RationalPair>;

RationalPair text(const cc::FT& value)
{
    auto t = cc::rational_to_text(value);
    return {std::move(t.numerator), std::move(t.denominator)};
}

PointText point_text(const cc::Point_2& p) { return {text(p.x()), text(p.y())}; }

cc::FT from_text(const RationalPair& value) { return cc::rational_from_text(value.first, value.second); }

std::vector<cc::VertexDisk> vertex_disks(const std::vector<std::vector<RationalPair>>& rows)
{
    std::vector<cc::VertexDisk> disks;
    disks.reserve(rows.size());
    for (const auto& row : rows) {
        if (row.size() != 5) throw cc::InvalidCertificateInputError("a vertex disk is (x0, x1, y0, y1, radicand)");
        disks.push_back(cc::VertexDisk::build(from_text(row[0]), from_text(row[1]), from_text(row[2]),
                                              from_text(row[3]), from_text(row[4])));
    }
    return disks;
}

cc::OneRootPoint one_root_point(const std::vector<RationalPair>& v)
{
    if (v.size() != 5) throw cc::InvalidCertificateInputError("a one-root point is (x0, x1, y0, y1, radicand)");
    return cc::OneRootPoint::build(from_text(v[0]), from_text(v[1]), from_text(v[2]), from_text(v[3]), from_text(v[4]));
}

cc::ExactTarget exact_target(const double tool_radius, const std::vector<PrimitiveRow>& primitive_rows,
                             const std::vector<PieceRow>& piece_rows, const std::vector<std::vector<RationalPair>>& disk_rows)
{
    std::vector<cc::DeclaredPrimitiveSpec> primitives;
    for (const auto& [is_line, sense, v] : primitive_rows) {
        if (v.size() != 7) throw cc::InvalidCertificateInputError("a declared primitive carries seven rationals");
        cc::DeclaredPrimitiveSpec p{is_line, cc::Point_2(from_text(v[0]), from_text(v[1])),
                                    cc::Point_2(from_text(v[2]), from_text(v[3])), cc::FT(0), cc::FT(0), cc::FT(0),
                                    cc::Point_2(0, 0), cc::FT(0), sense};
        if (is_line) {
            p.ux = from_text(v[4]);
            p.uy = from_text(v[5]);
            p.length = from_text(v[6]);
        } else {
            p.centre = cc::Point_2(from_text(v[4]), from_text(v[5]));
            p.radius = from_text(v[6]);
        }
        primitives.push_back(std::move(p));
    }
    std::vector<cc::ExactPieceSpec> pieces;
    for (const auto& [is_line, sense, support, source, target] : piece_rows) {
        cc::ExactPieceSpec p{is_line, cc::FT(0), cc::FT(0), cc::FT(0), cc::FT(0), cc::FT(0), cc::Point_2(0, 0), cc::FT(0),
                             cc::FT(0), sense, one_root_point(source), one_root_point(target)};
        if (is_line) {
            if (support.size() != 3) throw cc::InvalidCertificateInputError("an offset-line piece support is (nx, ny, level)");
            p.nx = from_text(support[0]);
            p.ny = from_text(support[1]);
            p.level = from_text(support[2]);
            p.ux = p.ny;
            p.uy = -p.nx;
        } else {
            if (support.size() != 4) throw cc::InvalidCertificateInputError("a circle piece support is (cx, cy, R, rho)");
            p.centre = cc::Point_2(from_text(support[0]), from_text(support[1]));
            p.radius = from_text(support[2]);
            p.boundary_radius = from_text(support[3]);
        }
        pieces.push_back(std::move(p));
    }
    return cc::ExactTarget::build(tool_radius, primitives, pieces, vertex_disks(disk_rows));
}

nb::dict counters_dict(const cc::DecisionCounters& c)
{
    nb::dict d;
    d["candidate_pairs"] = c.candidate_pairs;
    d["pairs_filtered_disjoint"] = c.pairs_filtered_disjoint;
    d["pairs_filtered_transversal"] = c.pairs_filtered_transversal;
    d["pairs_exact"] = c.pairs_exact;
    d["tangencies"] = c.tangencies;
    d["coincident_pairs"] = c.coincident_pairs;
    d["materialised_pairs"] = c.materialised_pairs;
    d["pair_towers"] = c.pair_towers;
    d["exact_branch_matches"] = c.exact_branch_matches;
    d["exact_half"] = c.exact_half;
    d["exact_compare_x"] = c.exact_compare_x;
    d["exact_membership"] = c.exact_membership;
    d["exact_target_membership"] = c.exact_target_membership;
    d["start_retries"] = c.start_retries;
    d["witness_bisections"] = c.witness_bisections;
    d["uncovered_halvings"] = c.uncovered_halvings;
    d["boundary_circles"] = c.boundary_circles;
    d["sub_arcs"] = c.sub_arcs;
    d["covered_by_coincident_support"] = c.covered_by_coincident_support;
    d["circle_line_pairs"] = c.circle_line_pairs;
    d["circle_line_exact"] = c.circle_line_exact;
    d["line_line_exact"] = c.line_line_exact;
    d["line_tangencies"] = c.line_tangencies;
    d["materialised_line_events"] = c.materialised_line_events;
    d["exact_tower_compares"] = c.exact_tower_compares;
    d["cross_base_separations"] = c.cross_base_separations;
    d["exact_line_side"] = c.exact_line_side;
    d["line_pieces"] = c.line_pieces;
    d["line_sub_edges"] = c.line_sub_edges;
    d["coincident_line_supports"] = c.coincident_line_supports;
    d["target_pieces"] = c.target_pieces;
    d["sweep_exposed_outside_target"] = c.sweep_exposed_outside_target;
    d["target_decisions"] = c.target_decisions;
    d["vertex_circle_pairs"] = c.vertex_circle_pairs;
    d["vertex_circle_exact"] = c.vertex_circle_exact;
    d["vertex_line_pairs"] = c.vertex_line_pairs;
    d["vertex_line_exact"] = c.vertex_line_exact;
    d["vertex_circle_tangencies"] = c.vertex_circle_tangencies;
    d["materialised_vertex_events"] = c.materialised_vertex_events;
    d["exact_target_events"] = c.exact_target_events;
    d["exact_fallbacks"] = c.exact_fallbacks();
    return d;
}
}  // namespace

NB_MODULE(_coverage_certificate_2, m)
{
    nb::exception<cc::InvalidCentreDomainError>(m, "InvalidCentreDomainError", PyExc_ValueError);
    nb::exception<cc::InvalidCertificateInputError>(m, "InvalidCertificateInputError", PyExc_ValueError);
    auto native_error = nb::exception<cc::CoverageCertificateNativeError>(m, "CoverageCertificateNativeError", PyExc_RuntimeError);
    // Subclasses after their base: nanobind tries the latest translator first.
    nb::exception<cc::StartWitnessSearchError>(m, "StartWitnessSearchError", native_error);
    nb::exception<cc::SignRuleContractError>(m, "SignRuleContractError", native_error);
    nb::exception<cc::SeparationBoundExceededError>(m, "SeparationBoundExceededError", native_error);
    nb::exception<cc::CoincidentTargetSupportError>(m, "CoincidentTargetSupportError", native_error);
    nb::exception<compas_cgal::exact::NonFiniteBinary64Error>(m, "NonFiniteBinary64Error", PyExc_ValueError);

    // Contract probes: the fail-loud bounds of the exact search loops, exercised
    // directly (tests/test_coverage_certificate_contract_bounds.py). Never part
    // of a certificate run.
    m.def(
        "contract_probe_start_witness",
        [](const double cx, const double cy, const double radius, const bool coincident) {
            const cc::FT r = cc::inject_binary64(radius);
            const cc::Point_2 c(cc::inject_binary64(cx), cc::inject_binary64(cy));
            const cc::Circle_2 circle(c, r * r);
            const cc::BoundaryCircle C{0, cc::PieceKind::Disk, circle, r, cc::interval_of(c.x()), cc::interval_of(c.y()),
                                       cc::interval_of(r * r), circle.bbox()};
            // One support: C itself (an uncut coincident support) or the circle of
            // radius r about c + (r, 0), which crosses C twice.
            const cc::Point_2 other = coincident ? c : cc::Point_2(c.x() + r, c.y());
            cc::DecisionCounters counters;
            const std::uint32_t bound = cc::start_witness_attempt_bound(0, 1, 0);
            const auto start = cc::find_start_witness(
                C, bound,
                [&](const cc::Point_2& w) { return CGAL::compare(CGAL::squared_distance(w, other), r * r) == CGAL::EQUAL; },
                counters);
            return std::make_pair(start.attempt + 1, bound);
        },
        "cx"_a, "cy"_a, "radius"_a, "coincident"_a);
    m.def(
        "contract_probe_sign_check",
        [](const int claimed, const double lo, const double hi) {
            (void)cc::require_sign_agrees(static_cast<CGAL::Sign>(claimed), cc::Interval(lo, hi), "probe");
        },
        "claimed"_a, "lo"_a, "hi"_a);
    m.def(
        "contract_probe_sign_rational_two",
        [](const RationalPair& u, const RationalPair& v, const RationalPair& p, const RationalPair& w, const RationalPair& q) {
            return static_cast<int>(cc::sign_rational_two(from_text(u), from_text(v), from_text(p), from_text(w), from_text(q)));
        },
        "u"_a, "v"_a, "p"_a, "w"_a, "q"_a);
    m.def(
        "contract_probe_rational_between",
        [](const RationalPair& lo_m, const RationalPair& lo_n, const RationalPair& lo_d, const RationalPair& hi_m,
           const RationalPair& hi_n, const RationalPair& hi_d) {
            const auto tower = [](const RationalPair& m_, const RationalPair& n_, const RationalPair& d_) {
                return cc::Tower{cc::OneRoot(from_text(m_)), cc::OneRoot(from_text(n_)), cc::OneRoot(from_text(d_))};
            };
            cc::DecisionCounters counters;
            const cc::ChartBound lo{false, cc::FT(0), tower(lo_m, lo_n, lo_d)};
            const cc::ChartBound hi{false, cc::FT(0), tower(hi_m, hi_n, hi_d)};
            return text(cc::rational_between(lo, hi, counters));
        },
        "lo_m"_a, "lo_n"_a, "lo_d"_a, "hi_m"_a, "hi_n"_a, "hi_d"_a);
    m.def(
        "contract_probe_vertex_disk",
        [](const RationalPair& x0, const RationalPair& x1, const RationalPair& y0, const RationalPair& y1,
           const RationalPair& delta) {
            const auto disk = cc::VertexDisk::build(from_text(x0), from_text(x1), from_text(y0), from_text(y1), from_text(delta));
            return std::make_pair(disk.x.is_extended(), disk.y.is_extended());
        },
        "x0"_a, "x1"_a, "y0"_a, "y1"_a, "delta"_a);

    m.def(
        "contract_probe_exact_target_contains",
        [](const double tool_radius, const std::vector<PrimitiveRow>& boundary, const std::vector<PieceRow>& pieces,
           const std::vector<std::vector<RationalPair>>& disks, const std::vector<PointText>& points, const bool exact_only) {
            const auto target = exact_target(tool_radius, boundary, pieces, disks);
            cc::DecisionCounters counters;
            std::vector<bool> inside;
            for (const auto& [x, y] : points) {
                const cc::Point_2 p(from_text(x), from_text(y));
                inside.push_back(target.contains(cc::ExactPoint::rational(p),
                                                 exact_only ? cc::FilterMode::ExactOnly : cc::FilterMode::Lazy, counters));
            }
            return inside;
        },
        "tool_radius"_a, "boundary"_a, "pieces"_a, "vertex_disks"_a, "points"_a, "exact_only"_a);

    m.def(
        "validate_centre_domain",
        [](const std::vector<XY>& outer, const std::vector<std::vector<XY>>& holes) {
            nb::gil_scoped_release release;
            (void)cc::CentreDomain::build(outer, holes);
        },
        "outer"_a, "holes"_a);

    m.def(
        "run_circle_core",
        [](const std::vector<XY>& centres, const std::vector<double>& path_radii, const double tool_radius,
           const std::vector<XY>& outer, const std::vector<std::vector<XY>>& holes,
           const std::vector<std::vector<RationalPair>>& disks, const bool exact_only, const bool audit) {
            cc::CircleCoreResult result;
            std::vector<EdgeRecord> records;
            {
                nb::gil_scoped_release release;
                const auto set = cc::CircleRegionSet::build(centres, path_radii, tool_radius);
                const auto target =
                    cc::MachinableTarget::build(cc::CentreDomain::build(outer, holes), tool_radius, vertex_disks(disks));
                result = cc::run_circle_core(set, target, exact_only ? cc::FilterMode::ExactOnly : cc::FilterMode::Lazy, audit);
                records.reserve(result.exposed.size());
                for (const auto& e : result.exposed) {
                    records.emplace_back(e.circle, e.region, static_cast<int>(e.kind), e.edge, e.whole_circle,
                                         e.start_vertex, e.end_vertex, text(e.parameter), point_text(e.on_edge),
                                         point_text(e.uncovered), e.uncovered_in_target);
                }
            }
            nb::dict out;
            out["exposed"] = records;
            out["counters"] = counters_dict(result.counters);
            out["region_count"] = result.region_count;
            out["piece_count"] = result.piece_count;
            out["elapsed_s"] = result.elapsed_s;
            out["phase_s"] = std::array<double, 4>{result.crossings_s, result.sweep_s, result.witness_s, result.target_s};
            return out;
        },
        "centres"_a, "path_radii"_a, "tool_radius"_a, "outer"_a, "holes"_a, "vertex_disks"_a, "exact_only"_a, "audit"_a);

    m.def(
        "run_certificate",
        [](const std::vector<XY>& centres, const std::vector<double>& path_radii, const double tool_radius,
           const std::vector<XY>& outer, const std::vector<std::vector<XY>>& holes,
           const std::vector<ConnectorRow>& connectors, const std::vector<std::vector<RationalPair>>& disks,
           const bool exact_only, const bool audit) {
            cc::CertificateResult result;
            std::vector<WitnessRecord> records;
            {
                nb::gil_scoped_release release;
                auto set = cc::CircleRegionSet::build(centres, path_radii, tool_radius);
                const auto target =
                    cc::MachinableTarget::build(cc::CentreDomain::build(outer, holes), tool_radius, vertex_disks(disks));
                std::vector<cc::ConnectorSpec> specs;
                specs.reserve(connectors.size());
                for (const auto& [loop, s0, t0, s1, t1] : connectors) specs.push_back({loop, s0, from_text(t0), s1, from_text(t1)});
                result = cc::run_certificate(std::move(set), target, specs,
                                             exact_only ? cc::FilterMode::ExactOnly : cc::FilterMode::Lazy, audit);
                records.reserve(result.exposed.size());
                for (const auto& e : result.exposed) {
                    records.emplace_back(static_cast<int>(e.kind), e.region_index, e.piece, e.edge, text(e.parameter),
                                         point_text(e.uncovered));
                }
            }
            nb::dict out;
            out["exposed"] = records;
            out["undecided"] = result.undecided;
            out["counters"] = counters_dict(result.counters);
            out["region_count"] = result.region_count;
            out["piece_count"] = result.piece_count;
            out["edge_count"] = result.edge_count;
            out["elapsed_s"] = result.elapsed_s;
            out["phase_s"] = std::array<double, 4>{result.crossings_s, result.circle_sweep_s, result.line_sweep_s,
                                                   result.witness_s};
            return out;
        },
        "centres"_a, "path_radii"_a, "tool_radius"_a, "outer"_a, "holes"_a, "connectors"_a, "vertex_disks"_a,
        "exact_only"_a, "audit"_a);

    m.def(
        "run_exact_target_certificate",
        [](const std::vector<XY>& centres, const std::vector<double>& path_radii, const double tool_radius,
           const std::vector<XY>& outer, const std::vector<std::vector<XY>>& holes,
           const std::vector<ConnectorRow>& connectors, const std::vector<PrimitiveRow>& boundary,
           const std::vector<PieceRow>& pieces, const std::vector<std::vector<RationalPair>>& disks, const bool exact_only,
           const bool audit) {
            cc::CertificateResult result;
            std::vector<WitnessRecord> records;
            {
                nb::gil_scoped_release release;
                auto set = cc::CircleRegionSet::build(centres, path_radii, tool_radius);
                const auto path_domain = cc::CentreDomain::build(outer, holes);
                const auto target = exact_target(tool_radius, boundary, pieces, disks);
                std::vector<cc::ConnectorSpec> specs;
                specs.reserve(connectors.size());
                for (const auto& [loop, s0, t0, s1, t1] : connectors) specs.push_back({loop, s0, from_text(t0), s1, from_text(t1)});
                result = cc::run_exact_target_certificate(std::move(set), path_domain, specs, target,
                                                          exact_only ? cc::FilterMode::ExactOnly : cc::FilterMode::Lazy, audit);
                records.reserve(result.exposed.size());
                for (const auto& e : result.exposed) {
                    records.emplace_back(static_cast<int>(e.kind), e.region_index, e.piece, e.edge, text(e.parameter),
                                         point_text(e.uncovered));
                }
            }
            nb::dict out;
            out["exposed"] = records;
            out["undecided"] = result.undecided;
            out["counters"] = counters_dict(result.counters);
            out["region_count"] = result.region_count;
            out["piece_count"] = result.piece_count;
            out["edge_count"] = result.edge_count;
            out["elapsed_s"] = result.elapsed_s;
            out["phase_s"] = std::array<double, 4>{result.crossings_s, result.circle_sweep_s, result.line_sweep_s,
                                                   result.witness_s};
            return out;
        },
        "centres"_a, "path_radii"_a, "tool_radius"_a, "outer"_a, "holes"_a, "connectors"_a, "boundary"_a, "pieces"_a,
        "vertex_disks"_a, "exact_only"_a, "audit"_a);
}
