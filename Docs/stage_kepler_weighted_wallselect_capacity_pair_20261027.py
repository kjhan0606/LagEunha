"""Reuse the accepted failure-point pair harness for the wall-select build."""

import stage_kepler_weighted_cache_pair_20261027 as pair

pair.BUILD=pair.ROOT/'gpu_energy_neutral_wallselect_cache_capacity64_build_v10'
pair.RUN=pair.ROOT/'gpu_energy_neutral_wallselect_capacity_pair_t006_v1'
pair.BINARY_NAME='eunha2_wallselect_cache_capacity64_v10.exe'

if __name__=='__main__':
    pair.main()
