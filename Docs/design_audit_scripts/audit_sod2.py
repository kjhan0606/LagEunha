import sys; sys.path.insert(0,'/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad')
from audit_run1d import sod
sod(hg=None, qlaw=1, zface=False, zh=False)             # cell Z, qlaw1, design dt
sod(hg=None, qlaw=1, zface=True, zh=False, dtmode='legacy')   # face Z only, legacy dt
sod(hg=None, qlaw=2, zface=False)                        # qlaw2 with cell Z
sod(dtmode='geom')                                       # design, geometric dt only
sod(dtmode=1.0)                                          # spectral coefficient 1.0
sod(dtmode=2.0)                                          # spectral coefficient 2.0
sod(hg=None, dtmode='geom')
