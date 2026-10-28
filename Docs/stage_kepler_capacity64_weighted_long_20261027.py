"""Continue the fresh weighted chain to t=10 and conditionally t=20."""

import stage_kepler_capacity64_weighted_chain_20261027 as chain

ROOT=chain.ROOT
P65=ROOT/'gpu_energy_neutral_capacity64_v5_recover_to65_a10_v1'
P10=ROOT/'gpu_energy_neutral_capacity64_v5_fresh_to10_a10_v1'
P20=ROOT/'gpu_energy_neutral_capacity64_v5_fresh_to20_a10_v1'
chain.RUNS=[(P10,P65,'audit_t65.json',6.5,10.,'audit_t10.json'),
            (P20,P10,'audit_t10.json',10.,20.,'audit_t20.json')]
chain.WALLTIME='06:00:00'

if __name__=='__main__':
    chain.main()
