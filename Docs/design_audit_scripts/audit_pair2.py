import sys; sys.path.insert(0,'/tmp/claude-10396/-home-kjhan-BACKUP-Eunha-A1-Claude-LagForce/84402936-5ffd-4ee7-a57e-1b2d07d92f79/scratchpad')
from audit_run1d import pair
w = sys.argv[1]
if w in ('1','8','64'):
    ratio = float(w)
    for M in (1, 10, 100, 1000):
        pair(M, ratio, T=(20.0 if M <= 10 else 4.0), dtmode='geom')
elif w == 'abl':
    for M in (10, 1000):
        T = 20.0 if M <= 10 else 4.0
        pair(M, 1.0, T=T, hg=None, dtmode='geom'); pair(M, 1.0, T=T, qon=False, dtmode='geom'); pair(M, 1.0, T=T, hg='acou', dtmode='geom'); pair(M, 1.0, T=T, zface=False, dtmode='geom')
elif w == 'design_dt':
    pair(10, 1.0, T=20.0); pair(1000, 1.0, T=4.0)
