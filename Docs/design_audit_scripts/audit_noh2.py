import sys; sys.path.insert(0,'/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad')
from audit_run1d import noh
noh(amp=0.1, hg=None, qlaw=1, zface=False, zmode='ZF', dtmode='legacy')   # reproduce old cell+Z_F baseline (0.035)
noh(amp=0.1, zmode='ZF')                                                   # design but Z_F amplitude on face states
noh(amp=0.1, zmode='ZH2')
noh(amp=0.1, zmode='davis2')
noh(amp=0.1, qlaw=1)
noh(amp=0.1, zface=False, zmode='ZF')                                      # cell-state Z_F with hg + qlaw2
noh(amp=0.1, hg=None, zface=False, zmode='ZF')                             # cell-state Z_F, qlaw2, no hg
noh(amp=0.1, hgcoef=4.0)
noh(dtmode='geom'); noh(amp=0.1, dtmode='geom')
