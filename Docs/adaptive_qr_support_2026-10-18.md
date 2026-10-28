# Continuous reconstruction support candidate

The fixed2.5sqrt(V) smooth support stopped before t4 because only8 samples
remained for9 cubic coefficients. Keep that failure and all old runs.

## Candidate definition

For distinct particle positions let T_j(R)=max(1-d_j^2/R^2,0)^2 and choose
R by sum_j T_j(R)=12. This is a kernel sum, not statistical effective sample
size. Use the existing LSQ row multiplier sqrt(V)/d_j times T_j(R), on both
design and RHS. These reconstruction weights are not Laguerre cell weights.

For positive distances, the sum is continuous and has derivative
4 sum_(d_j<R) d_j^2/R^3 (1-d_j^2/R^2)>0 at a positive target. With more
than12 distinct available points it increases from0 to a limit above12,
so its positive root is unique. The taper and its first derivative vanish
at the support boundary. Locally the implicit radius is continuously
differentiable while positions remain distinct and the full support is
available. No nearest-neighbor order switch is used in the radius equation.
Every nonzero term is below1, so the root requires more than12 samples.
This does NOT guarantee cubic rank for degenerate arrangements; the rank
guard remains necessary. The target12 is an experimental design choice,
not a proven optimal value or a physical constant.

## Independent saved-state audit

Full-cloud cKDTree audit, current r<1.9,11498 cells per state:

| State | Support count range | Maximum R/sqrt(V) | Maximum condition | Minimum rank |
|---|---:|---:|---:|---:|
| legacy t3.4004059 | 35–42 | 4.41287 | 6.91044 | 9 |
| legacy t4.0029059 | 33–44 | 4.91743 | 9.83211 | 9 |

The old fixed support has minimum counts11 and7, respectively. Kernel
sum residuals are below1.1e-14. Each query uses128 nearest non-self points
and verifies its root lies strictly inside the farthest queried distance;
thus omitted full-cloud points cannot contribute. Condition numbers use
the actual weighted cubic basis scaled by sqrt(V).

These are preserved legacy states, not the unavailable failed smooth state.
They establish feasibility on these geometries, not production support
coverage, force smoothness in evolution, boundary behavior, or stability.

## Implementation gate

Exam/exam.c currently obtains candidates from linked-cell neighborhoods,
expanding for Voronoi corner/security failures. That criterion alone does
not certify coverage of this larger reconstruction support. Before adopting
the radius, explicitly ensure the entire support disk lies inside searched
cells, expanding the candidate search as needed; unseen contributors must
not appear with nonzero weight when the cell neighborhood changes. Keep
the rank guard, validate exported polynomial moments and independent full
cloud support, and then check actual stepping and long physical evolution.
No production source change or replacement evolution was launched by this
audit. Future runs need intermediate restart outputs before the known
failure interval so this checkpoint gap does not recur.

Evidence: canonical_adaptive_qr_support_2026-10-18.json.
Analyzer: audit_adaptive_qr_support.py (run on lageunha, BLAS/OMP1).

## Opt-in implementation and build

Exam/exam.c now provides GFS_QR_ADAPTIVE=1, requiring GFS_QR_SMOOTH=1
and MPI1. It computes the guaranteed disk reach inside the clipped searched
linked-cell block, checks the kernel sum at that reach, and expands the
candidate block until the target root is covered. A55-iteration bisection
then obtains the radius; the existing QR uses that support/taper. Distinct
positions, candidate capacity and rank guards remain. No degree fallback
or force steering is added. This does not certify distributed halo coverage.

Separate build Codex_adaptive_qr_build_20261018 completed successfully on
lageunha. Original binaries and legacy evolution remain untouched. Initial
terminal probes are running under Codex_adaptive_qr_probes_20261018,
session59986: default-off plus, adaptive plus, adaptive minus, from existing
t3.4 states. These must be checked against independently computed full-cloud
support, gradients and cubic moments before any evolution is launched.
Scripts: run_adaptive_qr_initial_probes.py and
analyze_adaptive_qr_initial_probes.py. These probes do not replace the
historical physical failure interval or long-time validation.

Initial terminal probe outcome: default-off exited0; both adaptive probes
exited2. Both stop at owner36855, search radius14 linked cells,4290 candidates,
with QR_ADAPTIVE incomplete-support diagnostic. The clipped-grid disk gate
does not yet handle outer boundary coverage adequately and reaches the4096
capacity guard. This exposes the gap between the r<1.9 independent audit
and full-domain production use. No adaptive evolution was launched.
Next establish the physical/periodic/padded boundary contract and implement
a complete support certificate there. Do not simply increase capacity,
ignore missing coverage, or restrict the audit to interior cells and call
the full method validated. Existing failed probe directories are preserved.

## Periodic boundary correction and full-domain terminal validation

Source audit: KH uses paddingTreeVorork4Particles and mkLinkedList2D_oExam;
the padding routine supplies a finite-width periodic strip. Merely searching
more clipped grid cells cannot populate missing images beyond that strip.
The revised adaptive mode requires that periodic linked-list callback and
MPI1. It first tries the certified local disk. When insufficient, it directly
enumerates real particles using minimum-image periodic distances inside a
trial disk, expanding that disk until its kernel sum brackets12. The root
must remain below half the shorter box length; coincident points/capacity
violations abort. The resulting samples feed QR only; Voronoi cell/face
geometry is not rebuilt with the enlarged reconstruction stencil.

Separate build Codex_adaptive_qr_periodic_build_20261018 succeeded, SHA256
e033bcc1ffbd73df73bea2b916cf204b00d3b0f9a5987ccd9cd54f1970299e7c.
All three terminal probes under Codex_adaptive_qr_periodic_probes_20261018
exited0. Both +/- states were checked against independent periodic cKDTree
fits over ALL36864 cells, including outer/frozen cells:

- support mismatch count0, minimum cubic rank9;
- maximum independent unlimited-gradient error1.583e-13;
- maximum normalized cubic-moment error5.330e-15;
- limited/unlimited-gradient collinearity error<=5.778e-7 at rounded RHS
  precision,1711 limited velocity components per branch.

The raw exported-versus-RHS gradient difference37.15 is a limited/unlimited
comparison, not an operator accuracy failure (as in the prior smooth audit).
Default-off RHS differs from the older binary by at most2.711e-20, not
bitwise identical. Analyzer now explicitly rejects support mismatches.
These results validate this reconstruction on two saved states; actual
stepping, OMP/rejection behavior, conservation and long evolution remain
required. No replacement evolution has yet been launched. Direct fallback
cost and distributed-halo extension remain open implementation concerns.
Evidence: canonical_adaptive_qr_periodic_probes_2026-10-18.json.

## Actual steps, thread/rejection parity, and physical continuation

Four +/- adaptive steps at dt1e-4 and5e-5 exited0 with one finite accepted
step each and no retries. Inner half-squared departure response derivative
is2.96481465e-7 and2.96438422e-7, respectively: still positive. No matching
adaptive-generator closure or long-time stability follows from these values.
OMP1 versus OMP4, and forced rejected dt1e-4 followed by accepted5e-5 versus
direct5e-5, match all11 reader fields bitwise in both branches. Full binary
hashes differ. Forced-rejection branches each verify exactly one retry.
Raw adaptive hydro normalized momentum residual<=1.237e-15,total-energy
rate residual<=2.827e-17,interior thermodynamic identity relative L1
<=2.272e-13; no negative Q*du among221159 directed faces per state.
Scope remains raw hydro before gravity/HG/frozen-cell constraints.

Physical evolution launched on lageunha using unchanged preserved t2.0029059
reference checkpoint184 and the new adaptive binary; this is a changed
operator, not a restart of the failed fixed-support method. No evolved
fixed-support checkpoint existed to reuse. Existing legacy is preserved.
New root Codex_adaptive_qr_physical_to_t10_20261018/adaptive_from_t2;
runner session28624,timeout parent3267915,actual executable3267966.
MPI1/OMP4,dtcap.0025,target10.002905917763403. Restart dumps every.2 plus
final; estimated1012140976bytes with286011368144896bytes free at launch.
First observed4 accepted steps reach t2.0129059,finite1/retries0.
The run must cross the prior support failure at3.8 and physical failure
interval6-10; these short gates do not replace either requirement.

Evidence: canonical_adaptive_qr_steps_2026-10-18.json,
canonical_adaptive_qr_threads_2026-10-18.json,
canonical_adaptive_qr_retry_2026-10-18.json,
canonical_adaptive_qr_hydro_conservation_2026-10-18.json.
Runner: run_adaptive_qr_physical_to_t10.py.

## Historical failure interval: matched legacy t6 output available

Physical analyzer now accepts --root; comparator accepts --root and --case
adaptive_from_t2, retaining the existing legacy continuation for comparison.
Legacy checkpoint1783 is at actual time6.000405917763673. Inner background
(804 original-IC material cells) phase RMS11.47080deg,relative-pressure
RMS1.31701,radial RMS.0114211,circular-departure velocity RMS.0206796,
exact-position RMS.0758727,log-entropy RMS.547990. Its t4.0029059 values
were1.51835deg,.338565,.00359772,.00845082,.0103578,.0733208 respectively.
This is continuing physical degradation in the existing matched legacy,
not the older pre-constrained-HG baseline. Mass is unchanged,weights zero.
The adaptive run still has only its copied t2 source frame; at live step38
it reaches t2.0979059 with finite1/retries0. No adaptive improvement claim.
Evidence: canonical_legacy_t6_physical_2026-10-18.json.

## Matching adaptive generator check

Existing terminal forces/HG modes and completed single steps were analyzed
without rerunning evolution. Same-state hashes and initialization match.
Inner generator is2.961776927e-7; measured forward derivatives at1e-4/5e-5
are2.964814653e-7/2.964384218e-7. Relative discrepancies are+.102564% and
+.0880313%, not a simple halving. Shear-band differences similarly change
+.106354% to+.091511%. Central differences halve approximately
(-.090510% to-.043563%), but edge error grows from+.001553% to+.022807%.
Core signal is only~1e-15 and relative differences change sign.

Thus broad generator/finite-map closure is NOT established. The inner
absolute residual2.60729e-10 at smaller dt is smaller than the old hard-cutoff
discrepancy~1.95e-9, but that is not a proof of its cause or of stability.
Both computed and measured inner growth remain positive. Reversible pressure
and HG at the identical state are unchanged; adaptive hydro-minus-pressure
contribution is-9.87836505e-8, compared with the unchanged HG-4.37110299e-8.
Remaining non-smooth force components, geometry transitions and finite-map
effects are not separated by this comparison. Preserve this open discrepancy
while the actual physical run advances; do not report a full consistency pass.
Evidence: canonical_adaptive_qr_generator_2026-10-18.json,
canonical_adaptive_qr_generator_comparison_2026-10-18.json.
Analyzer: compare_adaptive_qr_generator.py.

## First evolved physical comparison: t2.2004059

Adaptive checkpoint263 was fully saved and analyzed at the same actual time
as legacy263. Both retain fixed masses and zero cell weights. Inner-background
phase RMS is.128528918deg versus legacy.128525046deg; relative-pressure RMS
.031065988 versus.031066209; departure velocity RMS.000480112447 versus
.000480029718. Adaptive/legacy ratios: phase1.00003012,pressure.99999289,
radial.99995711,velocity1.00017234,position.99999209,log-entropy.99979732.
There is no appreciable inner physical improvement at this early endpoint.
It is only~.1975 time after changing the operator on the inherited t2 state,
not evidence through the3.8 support failure or6-10 physical failure interval.
Disk-core entropy error is~1% lower but remains a tiny absolute signal;
other core metrics are nearly equal. Do not generalize that isolated ratio.
Evidence: canonical_adaptive_qr_t2p2_comparison_2026-10-18.json.

## t2.4004059 physical comparison

Adaptive checkpoint343 is complete and matches legacy actual time. Inner
phase RMS.165688779deg versus.165696337deg; relative-pressure RMS.038340997
versus.038343320; circular-departure velocity RMS.000628708769 versus
.000628270883. Adaptive/legacy ratios are.99995438,.99993942,1.00069697.
Both methods have increased physical error since t2.2; differences remain
small and mixed, with adaptive departure velocity~.0697% worse. No clear
inner stability improvement is established. Disk-core log-entropy error is
~1.78% lower, but that small absolute diagnostic does not offset the main
inner physical failure gate. Keep the running historical-time validation.
Evidence: canonical_adaptive_qr_t2p4_comparison_2026-10-18.json.

## t2.6004059 physical comparison

Checkpoint423 matches the legacy actual time. Inner adaptive/legacy error
ratios: phase.99988544,pressure.99984108,radial.99983524,departure velocity
1.00073702,position.99986948,log-entropy1.00032051. Adaptive absolute inner
phase RMS.211091deg,relative-pressure RMS.0485322,departure velocity
RMS.000872665. Small mixed changes persist; no meaningful stabilization is
demonstrated. In the initial .3-.5 shear band, velocity error is~.0975% higher
and log-entropy error~.1443% higher. These early results do not establish
behavior at the known3.8 support failure or6-10 physical failure interval.
Evidence: canonical_adaptive_qr_t2p6_comparison_2026-10-18.json.

## t2.8004059 physical comparison

Saved adaptive503 matches legacy actual time. Inner adaptive/legacy ratios:
phase.99982901,pressure.99977453,radial.99992875,departure velocity1.00046302,
position.99984506,log-entropy1.00117668. Adaptive phase RMS.268873deg,
relative-pressure RMS.0578922,departure velocity RMS.001251787. Physical
errors continue growing in both methods; differences remain small and mixed.
This is not evidence of stabilization. Keep the same existing trajectories
through the former support failure and historical physical failure interval.
Evidence: canonical_adaptive_qr_t2p8_comparison_2026-10-18.json.

## t3.0004059 physical comparison

Checkpoint583 matches legacy actual time. Inner adaptive/legacy ratios:
phase.99981214,pressure.99976806,radial.99999586,departure velocity1.00037458,
position.99985277,log-entropy1.00374011. Adaptive absolute phase RMS.347808deg,
relative-pressure RMS.0742821,departure velocity RMS.001702267. The errors
grow in both methods. Nearly unchanged orbital/pressure error and~.374%
higher log-entropy error do not establish stabilization; retain the planned
failure-interval and long-time tests. Source t2 already has inherited error.
Evidence: canonical_adaptive_qr_t3p0_comparison_2026-10-18.json.

## t3.2004059 physical comparison

Complete checkpoint663 matches the original legacy actual time. Inner
adaptive/legacy error ratios: phase .999716274, pressure .999963594,
radial .999890310, departure velocity 1.000152279, position .999778337,
log-entropy 1.008297882. Adaptive absolute phase RMS is .455208351 deg,
relative-pressure RMS .104904438, and departure velocity RMS .002204236.
In the initial .3-.5 shear band, departure velocity error is .0271% higher
and log-entropy error 1.0280% higher. Orbital errors remain nearly unchanged;
the increasing inner entropy difference warrants retaining this diagnostic.
This does not establish stabilization or identify the full failure mechanism.
Both compared snapshots retain fixed masses and zero cell weights.

At this inspection the actual adaptive process3267966 was live at
t3.2779059 (step510), and legacy process3186332 was live at t7.7754059
(step1509), both with finite accepted states and no retry on the latest step.
Continue these trajectories without restarting; adaptive must still cross
the former t3.813 support failure and historical t6-10 physics interval.
Evidence: canonical_adaptive_qr_t3p2_comparison_2026-10-18.json.

## t3.4004059 physical comparison

Complete checkpoint743 matches the legacy actual time; masses remain fixed
and cell weights zero. Inner adaptive/legacy error ratios: phase .999546551,
pressure .999789010, radial .999919505, departure velocity 1.000307118,
position .999599435, log-entropy 1.014168668. Adaptive phase RMS is
.596064800 deg, relative-pressure RMS .161003342, departure velocity RMS
.003161774, and log-entropy RMS .022257904. The initial .3-.5 shear band
has departure velocity error .0407% higher and log-entropy error 1.6210%
higher. The inner entropy difference has increased from .83% at t3.2 to
1.42%, while physical orbital errors continue to grow almost unchanged.
This candidate has not demonstrated stabilization. Preserve the running
trajectory through the earlier support failure and t6-10 physical gate;
these comparisons cannot establish the eventual outcome in advance.
Evidence: canonical_adaptive_qr_t3p4_comparison_2026-10-18.json.

## Legacy reaches t8: physical failure remains clear

Complete legacy checkpoint2591 has actual time8.001133936046068, with
unchanged masses and zero weights. Inner-background phase RMS is28.994121
deg (t6:11.470797), relative-pressure RMS3.016741 (t6:1.317010), departure
velocity RMS.026521895, position RMS.183432313, log-entropy RMS1.135315.
The initial .3-.5 shear band's phase RMS is32.908269deg and relative-pressure
RMS3.547788. This is substantial physical degradation despite finite accepted
states and small global energy drift. It does not certify the new candidate.

The legacy timestep is now below .0025, so this output time differs from
the earlier .0004059 snapshot phase. Future adaptive/legacy t8 comparisons
must retain actual timestamps; do not label nominal t8 snapshots exactly
matched or relax the comparator's time gate silently. Individual physical
errors against each snapshot's analytic orbit remain valid. Preserve the
existing trajectories and use an explicitly justified comparison if needed.
Evidence: canonical_legacy_t8_physical_2026-10-18.json.

## t3.6004059 physical comparison

Complete adaptive checkpoint823 matches legacy time. Inner adaptive/legacy
ratios: phase .999334566, pressure 1.000064151, radial 1.000023808,
departure velocity .997477082, position .999438758, log-entropy 1.015646865.
Adaptive phase RMS .785488136deg, relative-pressure RMS .242488414,
departure velocity RMS .004877748. Departure velocity improves by .2523%
relative to legacy, but both trajectories' absolute errors grow substantially
since t3.4; entropy error is 1.5647% higher. The shear-band entropy error is
1.6996% higher. These mixed differences do not demonstrate stabilization.
Masses remain fixed and weights zero. Still await the former t3.813 support
failure and historical t6-10 physical failure gates without restarting.
Evidence: canonical_adaptive_qr_t3p6_comparison_2026-10-18.json.

## t3.8004059 physical comparison

Complete adaptive checkpoint903 matches legacy time, with fixed masses and
zero cell weights. Inner adaptive/legacy ratios: phase .998996936, pressure
1.001265090, radial .999011109, departure velocity .993672285, position
.999038500, log-entropy 1.018061197. Absolute adaptive phase RMS1.066507210
deg, relative-pressure RMS .288829018, departure velocity RMS .007209398.
The .6328% lower departure velocity error coexists with .1265% higher
pressure error and 1.8061% higher entropy error. Absolute orbital errors
continue growing: this is not physical stabilization. This saved time is
still just before the former fixed-support stop at3.8129059; its existence
alone does not prove that failure time has been crossed.
Evidence: canonical_adaptive_qr_t3p8_comparison_2026-10-18.json.

## Adaptive trajectory passes the previous fixed-support stop time

Process3267966 was directly observed live at accepted step741, elapsed
1.8524999999999716 from source2.0029059177634028, hence actual time
3.8554059177633744. Latest accepted state is finite with retries0. This
passes the fixed-support smooth trajectory's last accepted3.8129059 and
subsequent insufficient-sample abort. The adaptive trajectory is different,
so this establishes continued execution across the old failure time, not
that every old geometry has adequate rank or that future rank loss is
impossible. Global mass drift is zero; energy drift from restart is
-5.2660140e-8. Physical orbital stability remains unproven, with the t3.8
comparison already showing substantial growing error. Continue to t6-10.
Legacy process3186332 is also live at actual8.4379159.
Evidence: canonical_adaptive_qr_crossed_old_stop_2026-10-18.log.

## t4.0004059 physical comparison

Complete adaptive checkpoint983 matches legacy983 at4.00040591776336
(not continuation source984 at4.0029059). Masses remain fixed, weights zero.
Inner adaptive/legacy ratios: phase .998356203, pressure 1.001317267,
radial .997313244, departure velocity .991416337, position .998191470,
log-entropy 1.024545510. Adaptive phase RMS1.509194622deg, pressure RMS
.338188473, departure velocity RMS .008360184, log-entropy RMS .074770400.
The .8584% velocity improvement remains small relative to the growing
absolute physical error; entropy error is2.4546% higher. The shear band has
velocity error .8770% lower and entropy error2.5047% higher. Reconstruction
survival through the previous stop is established, but orbital stabilization
is not. Continue the existing trajectory through t6-10. After t4, legacy
comparison snapshots are sparse; unmatched adaptive frames must remain
unmatched, while their individual analytic-orbit errors can still be audited.
Evidence: canonical_adaptive_qr_t4p0_comparison_2026-10-18.json.

## t4.2004059 individual physical audit

Complete adaptive checkpoint1063 retains fixed masses and zero weights.
Inner phase RMS2.110579033deg, relative-pressure RMS .429281078, departure
velocity RMS .010038020, position RMS .014395781, log-entropy RMS .104751422.
At t4 these were1.509194622deg, .338188473, .008360184, .010293089,
and .074770400. Thus all five errors increase over this interval. Shear-band
phase RMS2.353488827deg and relative-pressure RMS .456068510 also show
continued degradation. No same-time legacy snapshot is available here;
this is an original-IC analytic-orbit audit, not a cross-run ratio or
interpolated comparison. The candidate has not arrested physical error
growth by t4.2. Preserve the t6-10 trajectory and checkpoints.
Evidence: canonical_adaptive_qr_t4p2_physical_2026-10-18.json.

## t4.4004059 individual physical audit

Complete adaptive checkpoint1143 retains fixed masses and zero weights.
Inner phase RMS2.833682328deg, relative-pressure RMS .505327621, departure
velocity RMS .011801390, position RMS .019247503, log-entropy RMS .138851309.
All five increase from t4.2. The initial .3-.5 shear band has phase RMS
3.176450962deg and relative-pressure RMS .541312147. There is no matched
legacy frame at this time: these are analytic-orbit/initial-state errors,
not a cross-run comparison. Physical degradation continues after passing
the old reconstruction abort time; no stabilization is demonstrated.
Evidence: canonical_adaptive_qr_t4p4_physical_2026-10-18.json.

## Legacy retry observation near t9.079

Direct process observation found legacy3186332 live at accepted step2084,
actual time9.079039804586392, finite1 and retries1. The nearby retry log
reports bad=0, over=1: proposed dt1.8336e-3 exceeds updated limit1.8326e-3.
Thus this observed retry is a timestep-limit rejection, not evidence of a
nonfinite state. Subsequent nearby attempts also show small limit overshoots.
Do not continue describing the legacy run as universally retry-free.
Adaptive3267966 was live at accepted step993, actual4.485405917763361,
finite1 and retries0. Both existing trajectories are retained.

## t4.6004059 individual physical audit

Complete adaptive checkpoint1223 retains fixed masses and zero weights.
Inner phase RMS3.650051855deg, relative-pressure RMS .597511332, departure
velocity RMS .012825869, position RMS .024654242, log-entropy RMS .178214807.
All five increase from t4.4. Shear-band phase RMS4.078651647deg and
relative-pressure RMS .662785194 continue the localized degradation.
No same-time legacy frame exists here; these are individual analytic-orbit
and original-IC errors. Reconstruction survival has not arrested physical
error growth. Continue the preserved trajectory through the t6-10 gate.
Evidence: canonical_adaptive_qr_t4p6_physical_2026-10-18.json.

## t4.8004059 individual physical audit

Complete adaptive checkpoint1303 retains fixed masses and zero weights.
Inner phase RMS4.543508383deg, relative-pressure RMS .614433182, departure
velocity RMS .014603651, position RMS .030499477, log-entropy RMS .225418614.
All five increase from t4.6. The initial .3-.5 shear band reaches phase
RMS5.053952682deg and relative-pressure RMS .672452920. The smaller recent
pressure increment does not imply stabilization: orbital and entropy errors
continue growing. No same-time legacy frame is available. Continue the
existing trajectory to the historical t6-10 gate.
Evidence: canonical_adaptive_qr_t4p8_physical_2026-10-18.json.

## t5.0004059 individual physical audit

Complete adaptive checkpoint1383 retains fixed masses and zero weights.
Inner phase RMS5.469826236deg, relative-pressure RMS .743143618, departure
velocity RMS .014886298, position RMS .036563808, log-entropy RMS .276586697.
All five increase from t4.8. The initial .3-.5 shear band has phase RMS
6.052054866deg and relative-pressure RMS .815753220. No same-time legacy
frame exists at this time. The candidate has continued past the old support
abort, but has not stabilized the physical trajectory by t5. Preserve the
running t6-10 validation; this is not a long-time or backend consistency pass.
Evidence: canonical_adaptive_qr_t5p0_physical_2026-10-18.json.

## t5.2004059 individual physical audit

Complete adaptive checkpoint1463 retains fixed masses and zero weights.
Inner phase RMS6.464142249deg, relative-pressure RMS .886926904, departure
velocity RMS .015581531, position RMS .043071607, log-entropy RMS .323699588.
All five increase from t5. The initial .3-.5 shear band has phase RMS
7.126334010deg and relative-pressure RMS .980628592. No matched legacy
snapshot exists here. These original-IC physical errors continue to grow;
the reconstruction change has not stabilized the trajectory. Retain the
existing historical failure-interval validation and all checkpoints.
Evidence: canonical_adaptive_qr_t5p2_physical_2026-10-18.json.

## t5.4004059 individual physical audit

Complete adaptive checkpoint1543 retains fixed masses and zero weights.
Inner phase RMS7.572082684deg, relative-pressure RMS .939866940, departure
velocity RMS .017781591, position RMS .050305710, log-entropy RMS .379775651.
All five increase from t5.2. The initial .3-.5 shear band has phase RMS
8.353665246deg and relative-pressure RMS1.040393158. Thus pressure errors
in that band have reached order unity even though the integration continues.
No matched legacy snapshot is available here; these are original-IC errors.
The candidate has not stabilized the physical orbit. Preserve the existing
t6-10 validation and do not interpret finite-state survival as accuracy.
Evidence: canonical_adaptive_qr_t5p4_physical_2026-10-18.json.

## Legacy continuation completes t10 with physical failure

Actual process3186332 is absent and legacy_continue/exit.code is0. Final
accepted step2551 is finite with retries0 and elapsed6.0003986074967077;
complete final checkpoint3535 has actual time10.003304525260502. Masses
remain fixed and weights zero. Inner phase RMS55.415317911deg, relative
pressure RMS6.325156785, departure velocity RMS .038610002, position RMS
.318033458, log-entropy RMS1.714430941. Shear-band phase RMS63.380703398deg
and relative-pressure RMS7.395199273 show severe physical degradation.
Global energy drift from the t4 restart is only -1.70569325e-7, illustrating
that a small global conservation residual does not establish orbital accuracy.
The run completed execution but failed physical stability. Preserve it as
the completed legacy reference; do not rerun. Adaptive3267966 remains live.
Final actual times must be checked before any adaptive/legacy comparison:
the legacy endpoint overshoots its target and is not exactly nominal t10.
Evidence: canonical_legacy_t10_physical_2026-10-18.json; legacy log/exit.code.

## t5.6004059 individual physical audit

Complete adaptive checkpoint1623 retains fixed masses and zero weights.
Inner phase RMS8.778589058deg, relative-pressure RMS1.092573995, departure
velocity RMS .018225967, position RMS .058248758, log-entropy RMS .439747231.
All five increase from t5.4. Shear-band phase RMS9.728224764deg and pressure
RMS1.212685160 continue the degradation. Inner-background pressure error
is now also order unity. No same-time legacy frame exists here. Physical
stabilization is not established; retain the running historical t6-10 gate.
Evidence: canonical_adaptive_qr_t5p6_physical_2026-10-18.json.

## t5.8004059 individual physical audit

Complete adaptive checkpoint1703 retains fixed masses and zero weights.
Inner phase RMS10.082433830deg, relative-pressure RMS1.239594859, departure
velocity RMS .020142625, position RMS .066797344, log-entropy RMS .497233905.
All five increase from t5.6. Shear-band phase RMS11.216985033deg and relative
pressure RMS1.399901529 continue the physical degradation. No same-time
legacy snapshot exists here. The next nominal t6 frame can be compared
only after verifying its actual time against legacy6.000405917763673.
This result is not a physical stability pass.
Evidence: canonical_adaptive_qr_t5p8_physical_2026-10-18.json.

## t6.0004059 matched physical failure-interval comparison

Complete adaptive1783 and legacy1783 have identical actual time
6.000405917763673, unchanged masses, zero weights and matching material masks.
Inner phase RMS is11.423347905deg versus11.470796573 legacy (ratio .995863525).
Relative-pressure RMS1.374946918 versus1.317009568 (ratio1.043991594);
departure velocity RMS .021075650 versus .020679587 (ratio1.019152366);
log-entropy RMS .558544786 versus .547989969 (ratio1.019260967).
Radial and position ratios are .990363376 and .996124521. The shear-band
pressure error is4.8277% higher and velocity error1.9826% higher.

This historical failure-interval checkpoint shows no material orbital
stabilization: phase improvement is only .4136%, while pressure and velocity
errors are worse. Reconstruction continuity/support repair is not sufficient
to cure the observed physical failure. Keep the existing t10 continuation
for its authorized endpoint, but do not promote this candidate as a stable
Kepler method or use its execution survival as a backend-physics pass.
The unmatched t4.2-5.8 frames remain explicitly unmatched; no interpolation
or nominal-time substitution was used.
Evidence: canonical_adaptive_qr_t6p0_comparison_2026-10-18.json.

## t6.2004059 individual physical audit

Complete adaptive checkpoint1863 retains fixed masses and zero weights.
Inner phase RMS12.810891250deg, relative-pressure RMS1.509477816, departure
velocity RMS .020384287, position RMS .084598156, log-entropy RMS .616266906.
Departure velocity decreases from .021075650 at t6, whereas phase, pressure,
position and entropy errors increase. This mixed evolution does not establish
stabilization. It also illustrates why the positive instantaneous generator
at t6 cannot be extrapolated as a constant growth rate over a .2 interval.
No same-time legacy frame exists here. Direct observation finds process3267966
live at accepted step1690, actual6.2279059, finite1 and retries0.
Evidence: canonical_adaptive_qr_t6p2_physical_2026-10-18.json.
