# Preserved original baseline reaches t16: physical failure remains

Read the completed original baseline snapshot017503 at
t16.00003474705187, SHA256
7e4fcf9891d2029a305faac1fd87ea94b6f9865508f32a004fa21069b7022597.
The original process2286374 remains live towardt20; at observation
2026-10-17T21:43:14.970248+00:00 its accepted time was16.0053025.
Masses remain unchanged and weights zero in the saved state.

| Initial material region | Wrapped phase RMS(deg) | Relative pressure RMS | Radial RMS | Circular velocity RMS |
|---|---:|---:|---:|---:|
| Inner background |99.6623|44.8605|.0742206|.0825407|
| Inner disk edge |21.0347|46.3649|.0166841|.0151970|
| Disk core |.233775|.725359|.000633495|.000790465|

Pressure column is RMS(P/P_initial-1), not mean pressure or percentage.
Phase is wrapped modulo2pi, so large errors cannot count lost windings.
Exact-orbit position RMS is.480171 in the inner background, independently
confirming large orbit error. Positive/finite numerical states and small
total-energy drift are insufficient for physical accuracy.

This is the original baseline executable, predating constrained-HG fixes;
it is NOT the current matched legacy control and NOT the new smooth-QR
method. At the live observation its total energy drift was-3.81888e-7 and
angular-momentum drift-.00223029. The constrained-HG original process was
still live at15.7243283, and the from-t0 dt-capped process at13.0549936.
Their different-time results are not substituted for a matched t16 check.

The current smooth-vs-legacy physical experiment continues towardt10 from
preserved states. It must subsequently demonstrate physical accuracy over
longer times to meet the objective; this old baseline failure is not a
claim about how the new method will behave. Existing runs/checkpoints were
not restarted, deleted or replaced during this audit.

Evidence: canonical_original_baseline_t16_physical_2026-10-18.json.
Analyzer: analyze_smooth_qr_physical_snapshots.py --original-baseline-only.
