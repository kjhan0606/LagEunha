# Independent saved-state geometry at exact t=3

An independent all-site SciPy Voronoi construction, clipped to the same 1024-sided polygon at R=2.75, passed both auxiliary guard radii and the pre-existing 1e-6 relative-volume tolerance. Native versus independent cell volume differed by at most 5.9571e-8 across all 90,352 sites. This confirms the saved native geometry; it is not a continuum accuracy or boundary-convergence result.

The same original sites were also advanced analytically along their Plummer circular orbits, then independently tessellated. No native trajectory was modified, no hydrodynamic simulation was rerun, and no reference state was written as a restart.

| Initial region | Native mean absolute density change | Circular-reference geometric density change | Native density versus circular geometric reference |
| --- | ---: | ---: | ---: |
| Disk .5–2 | 0.0497882% | 0.0499734% | 0.0004171% |
| Core r<.1 | 4.95735% | 1.26566% | 4.66636% |
| Core .1–.5 | 0.203398% | 0.198998% | 0.0063950% |

All means are mass weighted and all regions use the original IDs/radii. The third column measures density changes of fixed-mass circular sites as their Voronoi volumes breathe. This is a useful kinematic reference, not an exact finite-cell Euler solution. The fourth column is a separate absolute mean, not an additive decomposition of the first two columns.

The disk's small density fluctuations closely follow kinematic mesh breathing at this time. The inner core has a much larger dynamical departure, consistent with the already identified timestep sensitivity but not proving its complete cause. The native t6.5/t10 checkpoints, common-restart timestep controls and boundary study remain required. This analysis took about 30 seconds on the Ada compute server and does not add work to each evolution step.
