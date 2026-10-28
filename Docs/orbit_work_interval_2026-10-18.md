# Continued early-growth work interval

Completed from the preserved instrumented checkpoint225 at
t=2.1054059177634006. Target is40 accepted steps of dt=.0025, ending near
2.2054059177634007. Existing trajectories have not been restarted.

Remote root: /gpfs/kjhan/LagForce/Codex_orbit_work_interval_20261018.
Runner session79130 terminated with exit0. All40 accepted steps were
finite and zero-retry; endpoint265 is t=2.2054059177633984.

Ran analyze_orbit_work_interval.py successfully on lageunha.
It verified all40 accepted steps with no retries,160 stage records and320
HG-bin records; checks matching IDs/masses and endpoint duration; reports
global closure and per-step signs. Its no-retry restriction deliberately
prevents rejected trial logs from contaminating the integral.

Source/config review:

- Input params set physical viscosity Re=0, Monaghan alpha=beta=0,
  and av_mode=5. Active log confirms av_mode5 and gpu_enabled0.
- TCF force uses normal pressure plus Q and a tangential pressure term.
  The usual non-TCF face force is gated by !tcf_ok.
- TCF rejects dynamic weights and schemeH in this current path.
- Split HG excludes the unsplit HG acceleration accumulation.
- Subzone is default-off and mutually exclusive with HG affine.
- Hyperviscosity is default0 in Params/header.c and absent from the input;
  no explicit runtime value has yet been exported by this diagnostic.

For now retain hydro_minus_TCF, rather than claiming this subtraction is
pure pressure globally. Non-TCF boundary faces and runtime defaults deserve
explicit accounting before narrowing the label. The present run provides
an operator budget, not a fix or a long-term physical validation pass.

## Completed result

Global departure energy increased from8.8421495075e-10 to9.9877128802e-10,
a factor1.12955711. Integrated contributions:

| Contribution | Increment |
|---|---:|
| Hydro minus TCF | +1.6982005093e-10 |
| Differential rotation | +6.7997309126e-11 |
| TCF | -8.6995970310e-11 |
| HG, finite measured increments | -3.6311813580e-11 |
| Actual net change | +1.1455633727e-10 |
| Closure residual | +4.6761110277e-14 |

Residual is0.04082% of net change. The signs above hold for every one of
the40 accepted steps. The measured direct viscous work is damping,
while hydro-minus-TCF and differential rotation supply departure energy.
This rules out direct net AV acceleration as the whole growth explanation
for this interval, but not indirect heating/pressure feedback. It also
does not demonstrate an unstable eigenmode: the prior velocity-pattern
response damped. Continued forcing and coupled position/thermodynamic
response remain possible.

Next priority is distinguishing entropy-generated pressure forces from
isentropic geometric pressure forces along the early trajectory, with
explicit pressure/non-pressure accounting. Preserve endpoint265 for
subsequent evolution. Half-dt interval convergence is not yet measured.
Evidence: canonical_orbit_work_interval_2026-10-18.json.
