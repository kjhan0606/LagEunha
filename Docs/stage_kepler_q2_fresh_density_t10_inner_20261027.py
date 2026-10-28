"""Stage the true fresh weighted t10 inner-ring density zoom."""

import stage_kepler_q2_fresh_density_t10_20261027 as base


base.DEST = base.ROOT / "fresh_weighted_vs_w0_density_t10_inner_cpu_v1"
base.NAMES = ("kepler_q2_fresh_weighted_vs_w0_density_t10_inner_cpu_20261027.slurm",
              *base.NAMES[1:])


if __name__ == "__main__":
    base.main()
