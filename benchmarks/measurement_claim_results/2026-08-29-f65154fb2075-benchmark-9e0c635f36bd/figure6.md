# Figure 6 reproduction — path length against engagement cap (rect_20x12)

Across the 4 of 8 caps where a constant-spacing path met the cap at all, the engagement-controlled path is on average **0.33x** the baseline's length — **shorter** — on rect_20x12 (tool 2). At the **2** cap(s) where BOTH paths met the cap it is **0.40x**, which is the fully like-for-like comparison in this table. The remaining **4** cap(s) have no compliant spacing at all: no trial spacing reached them, and neither did the controlled generator, whose lowest measured maximum over the sweep is 135.6 deg. The controlled path stayed under its own cap away from the entry at **2 of 8** caps, and the exact predicate never fired away from an entry at **3 of 8**. Where it did not, the ratio on that row compares a baseline selected for compliance against a controlled path that is not compliant, and the measured maximum in the row says by how much.

**controlled** is `engagement_controlled_toolpath`, asked for the cap on its row; its advance is regulated by the exact engagement predicate. **constant spacing** is `trochoidal_mat_toolpath_circular` at a fixed stepover, reproducing the MATHSM protocol: the generator has no cap input, so the cap enters only by selecting the shortest trial spacing that measured at or below it. Neither is one of the paper's curves. Both engagement figures are measured AFTER each chain's entry cut, which is a full slot for any generator entering solid stock and pins the raw maximum near a full turn regardless of the setting being swept. The exceedance column is an extra exact-predicate check carried for the controlled curve only; the baseline's compliance is already what selected it.

## Per-cap comparison

| cap (deg) | controlled length | controlled max TEA after entry (deg) | controlled meets cap | controlled exceedances after entry | spacing (tool diam.) | constant-spacing length | constant-spacing max TEA after entry (deg) | length ratio |
| ---: | ---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: |
| 20 | 8631.6 | 141.8 | no | 536 | — | no compliant spacing | — | — |
| 40 | 6539.2 | 141.8 | no | 152 | — | no compliant spacing | — | — |
| 60 | 3355.9 | 141.8 | no | 84 | — | no compliant spacing | — | — |
| 80 | 2242.9 | 141.8 | no | 44 | — | no compliant spacing | — | — |
| 100 | 1368.7 | 141.8 | no | 20 | 0.100 | 6823.6 | 98.7 | 0.20 |
| 120 | 1077.3 | 141.8 | no | 0 | 0.200 | 3497.1 | 117.8 | 0.31 |
| 140 | 930.9 | 135.6 | yes | 0 | 0.300 | 2198.6 | 138.4 | 0.42 |
| 160 | 844.7 | 154.6 | yes | 0 | 0.300 | 2198.6 | 138.4 | 0.38 |

## Constant-spacing trials

Spacing does not order engagement: the measured maximum falls and rises as the spacing widens, which is why the baseline is a minimum over every compliant trial rather than a bisection on spacing.

| spacing (tool diam.) | length | cut motions | max TEA after entry (deg) | raw max TEA (deg) |
| ---: | ---: | ---: | ---: | ---: |
| 0.025 | 24510.1 | 2673 | 131.1 | 360.0 |
| 0.050 | 12600.8 | 1345 | 110.1 | 360.0 |
| 0.075 | 7850.0 | 861 | 121.0 | 360.0 |
| 0.100 | 6823.6 | 697 | 98.7 | 360.0 |
| 0.125 | 4752.8 | 513 | 119.4 | 360.0 |
| 0.150 | 4072.7 | 441 | 131.5 | 360.0 |
| 0.200 | 3497.1 | 355 | 117.8 | 360.0 |
| 0.250 | 2471.7 | 259 | 128.8 | 360.0 |
| 0.300 | 2198.6 | 225 | 138.4 | 360.0 |
| 0.400 | 1718.4 | 171 | 166.5 | 360.0 |
| 0.500 | 1381.9 | 133 | 189.6 | 360.0 |
| 0.600 | 1243.1 | 115 | 205.7 | 360.0 |
