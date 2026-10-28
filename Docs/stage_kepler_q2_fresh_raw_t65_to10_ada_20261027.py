"""Continue the fully fresh raw q2 lineage from its accepted t=6.5 state."""

import stage_kepler_q2_raw_rate_t65_to10_ada_20261027 as base


base.SOURCE=base.ROOT/"gpu_energy_neutral_raw_rate_v11b_fresh_t1_to65_ada_v1"
base.DEST=base.ROOT/"gpu_energy_neutral_raw_rate_v11b_fresh_t65_to10_ada_v1"


if __name__=="__main__":
    base.main()
