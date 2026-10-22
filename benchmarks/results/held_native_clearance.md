Recorded 2026-10-22T17:13:07+00:00 at commit `42dcd289b2` (dirty tree) on Apple M1 Max.

**Differential audit.** Python oracle vs native oracle on the same path. A decision is
(element, index, verdict, reported primitive, primitive kind, segment); *listed* counts the elements that are not CLEAR.
Depths are compared bit for bit.

| workload | path | identical decisions | listed | depth bit mismatches | gouge circ / conn | touch circ / conn | Python s | native s | load (Python; native) |
|---|---|---|---:|---:|---|---|---:|---:|---|
| `figure5_refined` | polygon Ω | True | 1133 | 0 | 542 / 475 | 78 / 38 | 4.43 | 1.30 | 30–30; 30–29 |
| `figure5_refined` | Ω* | True | 1887 | 0 | 0 / 0 | 95 / 1792 | 5.56 | 0.62 | 32–32; 32–32 |
| `figure8_upper_lane` | polygon Ω | True | 1331 | 0 | 731 / 600 | 0 / 0 | 20.78 | 0.95 | 89–66; 66–66 |
| `figure8_upper_lane` | Ω* | True | 3978 | 0 | 0 / 0 | 0 / 3978 | 11.21 | 0.57 | 58–52; 52–52 |
| `figure8_upper_boundary` | polygon Ω | True | 491 | 0 | 254 / 237 | 0 / 0 | 8.16 | 0.43 | 49–47; 47–47 |
| `figure8_upper_boundary` | Ω* | True | 1270 | 0 | 0 / 0 | 0 / 1270 | 8.73 | 0.46 | 44–39; 39–39 |
| `figure8_crossed_skis_boundary` | polygon Ω | True | 1635 | 0 | 814 / 821 | 0 / 0 | 10.78 | 0.97 | 30–31; 31–31 |
| `figure8_crossed_skis_boundary` | Ω* | True | 1442 | 0 | 0 / 0 | 0 / 1442 | 5.64 | 0.39 | 30–28; 28–28 |
| `figure8_monstera_boundary_qualified_corners` | polygon Ω | True | 2941 | 0 | 1536 / 1364 | 20 / 21 | 11.19 | 1.28 | 41–37; 37–37 |
| `figure8_monstera_boundary_qualified_corners` | Ω* | True | 4013 | 0 | 0 / 0 | 31 / 3982 | 11.56 | 1.69 | 31–27; 27–25 |

The gouge and touch counts are the native oracle's; the identical-decisions column says the Python oracle's are the same.

**Generation A/B.** Ω* containment by each oracle in one process, 3 runs each, alternating order.
Byte-identical means one SHA-256 triple (initial path, refined circles, connectors) over every run.

| workload | byte-identical | circles | median generation s (Python / native) | Python / native | median refine s (Python / native) | load range |
|---|---|---:|---|---:|---|---|
| `figure5_refined` | True | 1804 | 4.70 / 11.21 | 0.42 | 0.83 / 0.77 | 24–284 |
| `figure8_crossed_skis_boundary` | True | 1442 | 6.55 / 6.86 | 0.95 | 1.36 / 0.35 | 74–163 |
| `figure8_monstera_boundary_qualified_corners` | True | 3983 | 17.20 / 15.29 | 1.12 | 3.12 / 1.26 | 14–68 |
