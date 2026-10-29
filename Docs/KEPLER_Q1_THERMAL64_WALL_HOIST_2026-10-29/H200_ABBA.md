# H200 ABBA qualification submitted

Job 439842 runs four sequential short t=10 to t=10.001 trials within one H200 allocation in baseline/candidate/candidate/baseline order. It uses the exact same accepted source checkpoint and equal timestep cap as the independent longer thermal64 comparison, with identical options apart from the executable. Four approximately 52 MB snapshots and small logs are expected. The job has no scheduler dependency; it was pending for GPU allocation quota at submission while job 439819 was running.

The prepared auditor requires exact non-pointer restart records, footer and numerical ledger parity, unchanged fixed masses and W0, positive tagged thermal payload, and zero reported energy injection. It measures median per-step wall time after two warmup steps and reports the observed ABBA ratio without declaring long-time physics, scaling or a required speedup. Four short samples are an initial timing measurement, not a production performance guarantee.
