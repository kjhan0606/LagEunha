# Validation within the existing allocation after t=3

The original 438957 CUDA+MPI+OpenMP evolution remains on its verified halo40
binary, with the same parameters, IC, target t=3 and native checkpoint
schedule. Its numerical input checksums and `stage.json` are unchanged.
The post-gate checker was extended, with its original source archived, to
call one finite auxiliary dispatcher after the original native MPI program
exits and both the native 15-check and radial 9-check gates pass at exact t=3.
It does not run during a native RK4 step.

This uses an existing reserved GPU after its original segment finishes.
The standalone seven-case request 439115 otherwise waits on the eight-GPU
user cap; syn103 also has only one physical A100 PCIe card, already occupied
by the original segment. No extra allocation, concurrent GPU simulation,
change to another project's jobs or login-node simulation is introduced.
Previously completed native one/two-rank Slurm tasks peaked at approximately
1.04/1.81 GiB per task, well below the parent's 32 GiB request. The dispatcher
requires its original two CPU threads, one node and at least 45 allocation
minutes remaining.

## Execution conditions and failure handling

The dispatcher checks its exact stage directory and parent job ID. It verifies
the current frozen batch, auditors, case wrappers and checksums. A stage-local
atomic lock prevents duplicate dispatch; a separate atomic benchmark claim
records which allocation executes the auxiliary work.

It proceeds only when 439115 is still unstarted and pending, no case outputs
exist, the workload directory and owner match, and no operator hold exists.
It holds and rechecks that request, creates its claim, sends a pending-state
cancellation filter to the Slurm controller, and confirms CANCELLED with zero
runtime. An already running/completed/partial standalone is left alone; an
ambiguous cancellation never starts an inline duplicate. If the check is
skipped, the original simulation continuation remains unchanged.

The frozen seven-case batch body then runs sequentially inside the parent's
allocation, using the same worker paths and MPI bootstrap as the standalone
batch. Each worker records the actual parent allocation and effective flags.
The physical checkpoint used by those short implementation comparisons
remains the preserved t=2.6549183918680868 state. Their purpose is operator
and implementation parity, not proof of long-time stability.

The auxiliary deadline is at most 40 minutes and is constrained by the actual
allocation end time with cleanup reserve. Process-group termination alone
is insufficient for Slurm-launched MPI descendants. The dispatcher also
queries numeric steps in this allocation after the original MPI exit and,
when needed, signals only those auxiliary steps. It never cancels the parent,
batch, extern or another job. These are the step-scoped operations documented
by [Slurm scancel](https://slurm.schedmd.com/scancel.html) and
[Slurm squeue](https://slurm.schedmd.com/squeue.html).
State rechecks during cleanup are 60 seconds apart. An actual read-only probe
confirmed `%i` step IDs on this Slurm 24.11 installation; `%T` is a job-only
format here and was removed before the hook ran.

A normal auxiliary failure or incomplete comparison is recorded separately
from the passing original physics checkpoint. After confirmed cleanup, the
original t=4 continuation proceeds. Unverified residual ranks or an operator
interruption cause an operational stop at the preserved checkpoint so that
the next run cannot overlap an unaccounted MPI process. The checkpoint is
not replayed by the dispatcher. All outcomes retain the original `next_job`.

## No automatic production promotion

Astra approved the inline-validation plan in principle and advised against
automatic promotion under the current performance contract. A 5% gain over
the new binary with options disabled does not establish a gain over the old
production executable, and equal static register/stack counts do not close
that gap. The exact single-rank combination, old-production timing and full
future restart/validation chain still require assessment from actual data.
The dispatcher therefore collects results and leaves the original binary,
physics and subsequent stages unchanged even if all seven cases pass.

Nine mocked tests cover gate/time/resource requirements, supported step-ID
format, no cancellation of an already running standalone, unconfirmed
cancellation, duplicate dispatch, normal auxiliary success/failure and refusal
to advance with unreaped steps. They make no real scheduler mutation or GPU
evolution and are not native runtime validation. The hook was installed while
the original simulation was still below t=3 and has not yet executed.

Authoritative artifacts:

- [Installation receipt](KEPLER_CHECKPOINT_VALIDATION_INSTALL_2026-10-29.json)
- [Slurm CLI revision and read-only probe](KEPLER_CHECKPOINT_VALIDATION_CLI_REVISION_2026-10-29.json)
- [Dispatcher source](run_kepler_checkpoint_validation_20261029.py)
- [Mocked contract tests](test_kepler_checkpoint_validation_20261029.py)
- [Installer](install_kepler_checkpoint_validation_20261029.py)

During execution the original `to3p0/post_checkpoint_validation.json` records
the phase and outcome. The benchmark directory's `inline_execution_claim.json`
identifies the actual executor. Its `qualifier.json`, `comparison.json`,
per-case native logs and dumps remain the numerical evidence. Parent Slurm
COMPLETED alone does not mean the auxiliary comparison passed. t>6 and the
GIZMO long endpoint remain unverified.
