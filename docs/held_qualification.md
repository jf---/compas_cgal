# Held qualification

!!! danger "Under repair: no published Held path is machining-qualified"
    The Held toolpaths in the [gallery](held_toolpath_gallery.md) and on the
    [performance](performance.md) and [Held performance](held_performance.md)
    pages are drafts. They are not fit to cut. Treat every length, circle count
    and generation time quoted for them as a draft measurement, not a result
    for a qualified path.

**Conclusion.** No Held path is machining-qualified: all five gouge. Coverage
is now certified both ways. Figure 5, crossed skis and upper boundary are
certified covered; Monstera is proven uncovered at two islands until its
certified repair, which is certified covered. The paths stay unqualified until
the exact centre domain Ω\* removes the gouges.

- **They gouge.** On all five published workloads, the tool disk leaves the
  analytic pocket. The deep gouges sit at concave arcs. This is decided
  exactly, path element by path element, in the
  [gouge audit](#exact-gouge-audit) below.
- **Monstera leaves uncut islands.** The generated Monstera path leaves two
  interior slivers of material uncut. The coverage certificate proves them with
  exact witnesses in two clusters. Select-certify-repair re-inserts the dropped
  source circles, and the repaired path is certified covered.

## Qualification overlays

Each overlay puts both verdicts on one generated path. The first caption line
is exact: the `CoverageCertificate.build` verdict on every circle and connector,
closing seam included, and the gouge-audit counts. Figure 5 is CERTIFIED
COVERED and gouges (deepest 0.021 mm, at a concave arc). Monstera is NOT
COVERED, with 8 exact witnesses at two slivers, A (45.48, 61.13) and
B (34.01, 71.42), and gouges (deepest 0.046 mm). The second caption line is
diagnostic: the swept fill, the uncut and gouged shading, the inset framing and
the inset areas are a raster at 0.05 mm per pixel. Insets frame the exact
witness clusters (ranked by witness count) and the deepest gouge.

![Figure 5 qualification overlay: certified covered, gouges at concave arcs](assets/images/held_figure5_qualification.png)

![Monstera qualification overlay: not covered at slivers A and B, gouges at concave arcs](assets/images/held_monstera_qualification.png)

Generator: `benchmarks/tools/held_qualification_overlay.py`. Regenerate with:

```bash
pixi run held-qualification-figures
```

The task generates both paths live (Monstera takes minutes), then runs the
certificate and the gouge audit before drawing.

## Coverage certificate

`CoverageCertificate.build` runs on the generated paths, closing seam connector
included. It decides coverage of T = Ω ⊕ D_r for the current
polygonal-projection Ω. A certified-covered row is therefore not a
machining qualification: the same paths gouge (next section). Method, theorem
and open review findings: [Held coverage certificate](held_coverage_certificate.md).

Generator: `tools/coverage_certificate_acceptance.py`. Regenerate with:

```bash
pixi run coverage-certificate-acceptance
```

--8<-- "benchmarks/results/coverage_certificate_acceptance.md"

## Why the paths gouge

The centre domain Ω the generators place circles in is an inward offset of the
polygonal projection of the pocket, not the exact erosion of the analytic
pocket. At a concave arc, each chord of the projection lies in material. A
circle tangent to such a chord therefore cuts into the wall by up to the chord
sagitta. [Held exact centre domain](held_exact_centre_domain.md#evidence-the-old-domain)
shows that the maximum gouge equals the maximum concave chord sagitta. The
gouges against straight lines are exact sub-ulp offsets of Ω's straight edges,
not geometric defects. The audit counts them and splits them out.

## Exact gouge audit

A path element gouges when `dist(element, ∂P) < r`. A touch, `dist == r`, is
allowed. The gouge/touch counts and the split by primitive are exact decisions.
Depths, and the split column above 10⁻⁹ mm, are reporting values in doubles.

Generator: `benchmarks/held_gouge_audit_table.py`. Regenerate with:

```bash
pixi run held-gouge-audit
```

The slow-tier test `benchmarks/acceptance/test_held_gouge_audit_figure5_slow.py`
pins the Figure 5 row. It records the defect: once Ω is the exact erosion, it
must fail and be replaced by `require_gouge_free` passing.

--8<-- "benchmarks/results/held_gouge_audit.md"

## What is certified and what is diagnostic

| claim | kind | evidence |
| --- | --- | --- |
| The listed elements gouge or touch | certified, exact predicate | gouge audit above |
| The Monstera path leaves uncut material | certified, exact witnesses | [coverage certificate](#coverage-certificate) above; slow test `benchmarks/acceptance/test_coverage_certificate_monstera_path_slow.py` |
| Figure 5, crossed skis, upper boundary and repaired Monstera cover T = Ω ⊕ D_r (polygonal Ω) | certified, `CoverageCertificate.build` returns a certificate | [coverage certificate](#coverage-certificate) above; repaired Monstera pinned by the same slow test |
| Coverage of the exact-Ω\* target | **not certified**; with `VERTEX_DISK` members the certificate raises `UnmodelledCertificatePieceError` instead of passing | [Held coverage certificate](held_coverage_certificate.md) |
| Gouge depths in mm | diagnostic, doubles | gouge audit above |
| Island areas and centroids | diagnostic, raster | [The Monstera evidence](held_coverage_certificate.md#the-monstera-evidence) |
| Engagement strips in the gallery | diagnostic, predecessor-circle model only | [gallery](held_toolpath_gallery.md) |
| Path lengths and generation times | draft measurements of unqualified paths | [Held performance](held_performance.md), [performance](performance.md) |

## The repair

Three designs close the gap. Each has a failing test to write before any code:

1. **Exact centre domain Ω\*.** It replaces the polygonal-projection offset with
   the exact erosion of the analytic pocket. This removes the concave-arc gouges.
   See [Held exact centre domain](held_exact_centre_domain.md).
2. **Coverage-preserving selection.** Boundary-order selection may drop a circle
   only if the certificate still holds. This restores the circles that covered
   the Monstera islands. See
   [Held coverage-preserving selection](held_coverage_preserving_selection.md).
3. **Coverage certificate, phase 2 (done).** Connector capsules and target
   pieces are in; coverage is certified, not only refuted, on the polygonal Ω.
   Still open: `VERTEX_DISK` pieces for the exact Ω\* and the two fail-loud
   fixes (C1, I1) in
   [Held coverage certificate](held_coverage_certificate.md).

The time budgets for these steps are on the [performance](performance.md) page.
A Held path counts as machining-qualified only when its gouge audit row has zero
gouging elements and `CoverageCertificate.build` returns a certificate for it.
