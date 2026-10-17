import json
from pathlib import Path
import continue_factor_gresho as driver
if __name__=='__main__':
    driver.SOURCES[8]=Path('/gpfs/kjhan/LagForce/Codex_D_finiteR_20261010_010542/corrected_loss/gresho_n8_dt0.001.checkpoint.json')
    for dt in (.001,.0005):print(json.dumps(driver.run(8,dt)),flush=True)
