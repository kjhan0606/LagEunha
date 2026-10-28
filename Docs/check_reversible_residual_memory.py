"""New auxiliary-energy candidate, linear 1D fixed lattice only.

No Kepler, nonlinear Euler, topology, shock, MPI or CUDA claim.
u_dot=v; m v_dot=m*c^2*D0^2*u-B.T*z; z_dot=B*v.
B=sqrt(m)*c/h*(I-neighbor_average), E_z=|z|^2/2.
"""
import json
import numpy as np

n=64; h=1/n; mass=h; c=1.
def d0(a): return (np.roll(a,-1)-np.roll(a,1))/(2*h)
def residual(a): return a-(np.roll(a,-1)+np.roll(a,1))/2
def B(a): return np.sqrt(mass)*c/h*residual(a)
def rhs(y):
    u,v,z=y
    return np.array([v,c*c*d0(d0(u))-B(z)/mass,B(v)])
def energy(y):
    u,v,z=y
    return .5*(mass*(v@v+c*c*(d0(u)@d0(u)))+z@z)
rng=np.random.default_rng(20261018)
y=rng.normal(size=(3,n));dy=rhs(y)
terms=[mass*y[1]@dy[1],mass*c*c*d0(y[0])@d0(dy[0]),y[2]@dy[2]]
closure=float(sum(terms)/sum(abs(t) for t in terms))
assert abs(closure)<1e-13
rows=[]
for k in [1,8,32]:
    q=2*np.pi*k/n
    omega_old=c/h*abs(np.sin(q))
    omega=np.sqrt(omega_old**2+(c/h*(1-np.cos(q)))**2)
    assert abs(omega-2*c/h*abs(np.sin(q/2)))<1e-12
    mode=np.cos(q*np.arange(n));mode/=np.linalg.norm(mode)
    # Start from a pure velocity seed: u=z=0, so no energy is invented.
    initial=np.array([np.zeros(n),mode,np.zeros(n)])
    exact_v=mode*np.cos(omega)
    runs=[]
    for dt in [.001,.0005]:
        y=initial.copy();steps=round(1/dt)
        for _ in range(steps):
            a=rhs(y);b=rhs(y+.5*dt*a);cc=rhs(y+.5*dt*b);d=rhs(y+dt*cc)
            y+=dt/6*(a+2*b+2*cc+d)
        runs.append({'dt':dt,'relative_energy_error':float(energy(y)/energy(initial)-1),
                     'velocity_L2_error':float(np.linalg.norm(y[1]-exact_v)),
                     'auxiliary_energy_fraction':float(.5*(y[2]@y[2])/energy(initial))})
    rows.append({'mode':k,'old_frequency':float(omega_old),'new_frequency':float(omega),'runs':runs})
scaling=[]
for nn in [32,64,128,256]:
    hh=1/nn;q=2*np.pi/nn
    scaling.append({'n':nn,'fixed_physical_mode_added_frequency_squared':float((c/hh*(1-np.cos(q)))**2)})
print(json.dumps({'scope':__doc__,'relative_instantaneous_energy_closure':closure,
                  'modes':rows,'smooth_mode_scaling':scaling},indent=2))
