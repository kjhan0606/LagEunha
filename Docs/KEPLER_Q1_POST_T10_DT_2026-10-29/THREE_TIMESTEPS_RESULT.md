# Three same-parent timestep controls at t=10 to 10.1

Job 439375 completed with exit status 0. All three runs use the exact same accepted t=10 checkpoint, frozen one-rank H200 executable, fixed mass and W0. Each reached t=10.1 and passed its 28 original checkpoint checks; actual step sizes are separated as required. These checks are not a core-accuracy criterion or a long-time stability proof. Full details and endpoint hashes are in `THREE_TIMESTEPS_COMPARISON.json`.

| Control | Accepted steps | Median timestep | Core r0<0.1 mean absolute incremental log entropy | Disk 0.5<=r0<2 mean absolute incremental log entropy |
|---|---:|---:|---:|---:|
| Coarse | 258 | 3.97564e-4 | 2.95535e-4 | 3.01252e-6 |
| Half cap | 2101 | 4.76003e-5 | 4.23066e-6 | 6.02580e-5 |
| Quarter cap | 4202 | 2.38001e-5 | 5.44525e-6 | 1.07603e-4 |

The naming refers to caps relative to the parent's recent minimum step, not equal halving of the actual coarse median. Half/quarter median ratio is exactly 0.5; half/coarse median ratio is 0.11973. The core improves strongly from coarse to half, then worsens modestly at quarter. The disk absolute incremental entropy increases as steps are reduced: quarter is about 1.79 times half and about 35.7 times coarse. The signed disk means are -7.67e-7, -6.24e-6 and -8.66e-6 respectively; the absolute measure must not be described as net physical heating.

This is a short comparison starting from an already damaged state. Smaller steps increase the number of float updates, but these data alone do not distinguish accumulated rounding from evolving geometry, real pressure work, or other numerical effects, and they do not establish an asymptotic order. A separate same-parent H200 double-thermal comparison is running under job 439819 to test that attribution; it is not an inherited-error repair.
