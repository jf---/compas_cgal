# Benchmark corpus result

**42 of 42** measured instances put the tool over its engagement cap, at **3372** cut motions in total. A further **0** could not be proved under the cap without any exceedance being demonstrated. The slowest instance is **islands_4x4** at **1217.33 s** to certify (64658 stations); across all 42 instances certification costs **1379x** generation (4272.6 s against 3.10 s). 0 instances failed to measure.

**truly exceeding** counts cut motions where the exact engagement predicate fired at a sampled cutter position: a demonstrated LOWER BOUND on how many motions are over the cap, never a certificate. **uncertified** counts operations whose cap could not be *proved*, which includes every operation the certifier declined to measure because its guard would not close; it is sound and over-counts. Neither number estimates the other.

## Per-instance measurements

| instance | family | generate (s) | certify (s) | stations | max TEA (deg) | truly exceeding | uncertified | unresolved | arr. vertices | max digits |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| islands_4x4 | topology | 0.095 | 1217.33 | 64658 | 360.0 | 394 | 2616 | 2203 | 26977 | 129 |
| islands_3x3 | topology | 0.073 | 654.12 | 44486 | 360.0 | 323 | 1788 | 1441 | 20820 | 129 |
| disk_r8 | analytic | 0.396 | 466.78 | 99328 | 253.4 | 37 | 4736 | 4699 | 14561 | 261 |
| islands_2x2 | topology | 0.034 | 327.00 | 25132 | 360.0 | 304 | 968 | 650 | 17021 | 129 |
| ngon_k128 | complexity | 0.353 | 276.34 | 72448 | 240.3 | 27 | 3329 | 3302 | 10700 | 255 |
| arcfrac_n16_r1 | complexity | 0.503 | 263.48 | 61926 | 360.0 | 49 | 2805 | 2756 | 7354 | 133 |
| islands_1x1 | topology | 0.017 | 205.36 | 13650 | 360.0 | 237 | 506 | 266 | 15933 | 129 |
| arcfrac_n16_r0.5 | complexity | 0.170 | 122.05 | 37241 | 360.0 | 52 | 1694 | 1641 | 6560 | 248 |
| stadium_l20_w3 | analytic | 0.319 | 95.45 | 37412 | 360.0 | 51 | 1615 | 1564 | 7531 | 259 |
| arcfrac_n16_r0.25 | complexity | 0.070 | 61.31 | 25729 | 360.0 | 43 | 1186 | 1143 | 5994 | 247 |
| arcchan_r10_w2_s90 | analytic | 0.436 | 56.76 | 35151 | 218.7 | 72 | 1643 | 1571 | 6080 | 244 |
| dumbbell_p4 | necks | 0.016 | 30.36 | 9170 | 360.0 | 58 | 388 | 330 | 4733 | 135 |
| dumbbell_p1.05 | necks | 0.016 | 27.46 | 9402 | 200.0 | 98 | 362 | 264 | 5153 | 133 |
| dumbbell_p1.1 | necks | 0.016 | 26.68 | 9328 | 200.1 | 89 | 361 | 272 | 5152 | 135 |
| dumbbell_p2.5 | necks | 0.017 | 26.31 | 9274 | 360.0 | 66 | 392 | 326 | 4959 | 133 |
| dumbbell_p3 | necks | 0.017 | 25.88 | 9086 | 360.0 | 62 | 382 | 320 | 4826 | 131 |
| half_turn_arm | degeneracy | 0.011 | 25.75 | 6709 | 360.0 | 43 | 261 | 214 | 4834 | 137 |
| pinch_exactly_tool | degeneracy | 0.019 | 25.03 | 8838 | 347.0 | 94 | 354 | 254 | 5922 | 135 |
| dumbbell_p2 | necks | 0.017 | 24.76 | 9310 | 352.8 | 67 | 394 | 327 | 5030 | 131 |
| dumbbell_p1.7 | necks | 0.017 | 22.82 | 9168 | 268.5 | 69 | 388 | 319 | 5077 | 135 |
| dumbbell_p1.4 | necks | 0.017 | 21.75 | 9126 | 226.9 | 83 | 371 | 288 | 5106 | 131 |
| tangent_island | degeneracy | 0.014 | 20.45 | 6088 | 360.0 | 89 | 235 | 145 | 5107 | 129 |
| dumbbell_p1.2 | necks | 0.028 | 20.31 | 8974 | 208.7 | 81 | 357 | 276 | 5051 | 135 |
| collinear_run | degeneracy | 0.008 | 19.44 | 5189 | 191.8 | 68 | 205 | 134 | 4926 | 129 |
| arcfrac_n16_r0 | complexity | 0.037 | 17.99 | 15776 | 249.1 | 47 | 752 | 705 | 5510 | 251 |
| ngon_k32 | complexity | 0.061 | 16.29 | 18112 | 240.0 | 27 | 864 | 837 | 3456 | 261 |
| cocircular_square | degeneracy | 0.075 | 16.01 | 4164 | 191.8 | 63 | 168 | 102 | 3904 | 129 |
| prec_k12_r10_d0 | precision | 0.021 | 12.59 | 12112 | 243.2 | 48 | 575 | 527 | 4988 | 138 |
| rect_20x10 | analytic | 0.006 | 12.20 | 4629 | 360.0 | 32 | 176 | 144 | 4422 | 139 |
| prec_k12_r1_d6 | precision | 0.019 | 12.11 | 11657 | 238.4 | 47 | 552 | 505 | 4888 | 145 |
| prec_k12_r0.01_d6 | precision | 0.018 | 11.98 | 11657 | 238.4 | 47 | 552 | 505 | 4898 | 152 |
| prec_k12_r10000_d6 | precision | 0.019 | 11.95 | 11657 | 238.4 | 47 | 552 | 505 | 4912 | 127 |
| prec_k12_r100_d6 | precision | 0.019 | 11.91 | 11657 | 238.4 | 47 | 552 | 505 | 4872 | 136 |
| prec_k12_r10_d12 | precision | 0.019 | 11.91 | 11657 | 238.4 | 47 | 552 | 505 | 4934 | 139 |
| prec_k12_r10_d6 | precision | 0.019 | 11.85 | 11657 | 238.4 | 47 | 552 | 505 | 4928 | 140 |
| prec_k12_r10_d9 | precision | 0.019 | 11.83 | 11657 | 238.4 | 47 | 552 | 505 | 4876 | 140 |
| prec_k12_r10_d4 | precision | 0.018 | 11.80 | 11657 | 238.4 | 47 | 552 | 505 | 4878 | 139 |
| prec_k12_r10_d2 | precision | 0.019 | 11.76 | 11657 | 238.5 | 47 | 552 | 505 | 4896 | 140 |
| prec_k12_r10_d1 | precision | 0.019 | 11.76 | 11640 | 238.9 | 49 | 551 | 502 | 4906 | 141 |
| ngon_k3 | complexity | 0.004 | 8.86 | 3513 | 202.4 | 66 | 120 | 54 | 3055 | 252 |
| ngon_k12 | complexity | 0.019 | 4.27 | 6792 | 238.0 | 27 | 324 | 297 | 2391 | 266 |
| ngon_k6 | complexity | 0.008 | 2.60 | 4020 | 219.9 | 34 | 174 | 140 | 2193 | 253 |
