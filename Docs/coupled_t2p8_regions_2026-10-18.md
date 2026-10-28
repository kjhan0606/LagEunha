# Continued response localizes in the inner shear band

Analyzed completed common checkpoint503 at t2.8004059177633858, duration
.7975 from the t2 seed. Circular-departure squared response ratio is8.35267
in hydro versus.979118 in gravity only, a quotient8.53081. Gravity dt
refinement changes its ratios by at most3.44e-10; hydro even/odd velocity
ratio is.003403. These remain finite-time diagnostic norms, not a proven
unstable eigenmode or a converged long-time result.

New fixed-initial-radius material-bin analysis of the hydro response:

| Initial radius | Departure-response squared growth | Share of final global departure response |
|---|---:|---:|
| 0-.1 | .183 | .033% |
| .1-.2 | .504 | .902% |
| .2-.3 | 3.932 | 7.421% |
| .3-.4 | 15.486 | 45.001% |
| .4-.5 | 8.964 | 46.250% |
| .5-.6 | 7.698 | .393% |

Thus91.25% lies in .3-.5, while the innermost velocity response damps.
Top10 cells account for17.76% of departure response and23.54% of acoustic
response; it is not entirely one exceptional pair. Acoustic response is
also concentrated in .3-.5 (about71.47%), but grows even in some inner
bins where velocity response decreases.

Existing initial radial-noise audit finds shear |r*dOmega/dr| RMS about
.715,2.494,4.534,5.421,4.879 over successive .1 bins through r=.5.
This spatial association is consistent with the strong-shear band in the
softened potential; it does not isolate shear as the sole cause. In
particular initial geometric Vdot/V RMS peaks in .2-.3, while the later
response is concentrated farther out. Seed generation and amplification
need not peak at the same radius. The two audits use their own documented
radius calculations and should not be treated as identical cell lists at
floating-point bin boundaries.

This qualifies the shorthand 'noise grows inward': current evidence does
not show monotonic growth toward the center. It concerns this prescribed
coupled seed, not a complete spectrum or the unperturbed trajectory's
entire error budget. Refined branches are still running toward t2.6004059;
do not transfer shorter-horizon convergence statements to this result.

Evidence: canonical_coupled_t2p8_comparison_2026-10-18.json,
canonical_coupled_t2p8_regions_2026-10-18.json,
canonical_initial_radial_noise_2026-10-18.json.
Script: audit_coupled_response_regions.py --step503.
Maintain sessions62605 and93108. Long-time physics remains unresolved.
