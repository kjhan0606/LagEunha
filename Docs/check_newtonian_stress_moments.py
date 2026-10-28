"""Small analytic fixtures, executed on lageunha, not the login node."""
import json
import numpy as np
# Module source is supplied ahead of this script over SSH stdin.
assert viscosity(KEPLER_REYNOLDS,1.,1.)==0.
polygon=np.array([[.2,-.3],[1.7,.1],[.4,1.4]])
V=.5*np.sum(polygon[:,0]*np.roll(polygon[:,1],-1)-polygon[:,1]*np.roll(polygon[:,0],-1))
rows=[]
for name,B,nu in [('rigid',np.array([[0.,-2.],[2.,0.]]),.3),
                  ('shear',np.array([[0.,2.],[0.,0.]]),.3),
                  ('kepler_infinite_Re',np.array([[.1,2.],[-.7,-.1]]),viscosity(KEPLER_REYNOLDS,1.,1.))]:
    tau=stress(B,1.,nu);total=volume_rates(B,tau,V)
    for a,b in zip(polygon,np.roll(polygon,-1,axis=0)):
        edge=b-a;length=np.linalg.norm(edge);n=np.array([edge[1],-edge[0]])/length
        x=(a+b)/2;v=np.array([.4,-.8])+B@x
        rate=face_rates(x,v,n,tau,length)
        opposite=face_rates(x,v,-n,tau,length)
        for k in total:
            assert np.array_equal(rate[k],-opposite[k])
            total[k]+=rate[k]
    heating=total['E']-.5*np.trace(total['Pi'])
    expected=V*np.sum(tau*B)
    defects={k:float(np.max(np.abs(total[k]))) for k in ('p','D','Pi')}
    assert max(defects.values())<2e-14
    assert abs(heating-expected)<2e-14
    if name!='shear':assert np.array_equal(tau,np.zeros((2,2))) and heating==0.
    rows.append(dict(case=name,defects=defects,heating=heating,expected=expected))
print(json.dumps(dict(passed=True,results=rows,scope='Constitutive law and affine constant-stress quadrature only; no evolution or production backend validation.'),indent=2))
