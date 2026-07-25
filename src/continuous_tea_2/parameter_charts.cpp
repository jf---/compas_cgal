#include "parameter_charts.h"

#include "../exact_algebraic_1.h"

#include <algorithm>
#include <cstdint>
#include <string_view>
#include <utility>

namespace {

std::string u64_record(std::size_t value)
{
    std::string result;
    for (int shift = 56; shift >= 0; shift -= 8) {
        result.push_back(
            static_cast<char>(
                (static_cast<std::uint64_t>(value) >> shift)
                & 0xffU));
    }
    return result;
}

std::string ccan_node(char kind, const std::string& payload)
{
    std::string result("CCAN\0\1", 6);
    result.push_back(kind);
    result += u64_record(payload.size());
    result += payload;
    return result;
}

std::string ccan_bytes(const std::string& value)
{
    return ccan_node('B', value);
}

std::string ccan_tagged(
    std::string_view tag,
    const std::string& payload)
{
    return ccan_node(
        'T',
        ccan_bytes(std::string(tag))
            + ccan_bytes(payload));
}

std::string seam_id(
    std::string_view family,
    std::size_t ordinal)
{
    return ccan_tagged(
        "parameter-chart-seam-v1",
        ccan_tagged(
            family,
            std::to_string(ordinal)));
}

bool charts_equal(
    const ParameterChart2& left,
    const ParameterChart2& right)
{
    return left.chart_id == right.chart_id
        && left.family == right.family
        && left.domain_low == right.domain_low
        && left.domain_high == right.domain_high
        && left.orientation == right.orientation
        && left.start_seam_id == right.start_seam_id
        && left.end_seam_id == right.end_seam_id
        && left.owns_start_seam == right.owns_start_seam
        && left.owns_end_seam == right.owns_end_seam;
}

} // namespace

std::vector<ParameterChart2> parameter_charts()
{
    std::vector<ParameterChart2> charts;
    charts.reserve(6);
    for (std::size_t index = 0; index < 4; ++index) {
        charts.push_back(
            {
                "center-quarter-" + std::to_string(index)
                    + "-v1",
                "center-circle",
                "0",
                "1",
                "ccw",
                seam_id("center-circle", index),
                seam_id("center-circle", (index + 1) % 4),
                true,
                false,
            });
    }
    for (std::size_t index = 0; index < 2; ++index) {
        charts.push_back(
            {
                "rim-half-" + std::to_string(index) + "-v1",
                "cutter-rim",
                "-1",
                "1",
                "ccw",
                seam_id("cutter-rim", index),
                seam_id("cutter-rim", (index + 1) % 2),
                true,
                false,
            });
    }
    return charts;
}

VerifiedEventPartition2 verify_chart_coverage(
    const std::vector<ParameterChart2>& charts)
{
    const std::vector<ParameterChart2> expected =
        parameter_charts();
    const bool valid =
        charts.size() == expected.size()
        && std::equal(
            charts.begin(),
            charts.end(),
            expected.begin(),
            charts_equal);

    EventPartitionCertificate2 certificate;
    certificate.build_evidence =
        exact_algebraic_backend_evidence();
    certificate.charts = charts;
    certificate.source_kind = "parameter-chart-coverage-v1";
    if (valid) {
        certificate.seams.reserve(charts.size());
        for (const ParameterChart2& chart : charts) {
            certificate.seams.push_back(
                {
                    chart.start_seam_id,
                    chart.chart_id,
                });
        }
    }
    return {
        valid
            ? ContinuousTeaVerdict::CERTIFIED
            : ContinuousTeaVerdict::UNRESOLVED_DEGENERACY,
        std::move(certificate),
    };
}
