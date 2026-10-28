"""Algebraic HG boundary-transfer identity; not a solver/stability test."""
import math
import random
import unittest

class BoundaryBalance(unittest.TestCase):
    def test_free_and_fixed_thermal_owner(self):
        rng = random.Random(20261017)
        for fixed_owner in (False, True):
            for case in range(500):
                masses = [rng.uniform(.5, 2) for _ in range(4)]
                modes = [(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in masses]
                modes[-1] = tuple(-sum(p[k] for p in modes[:-1]) for k in (0, 1))
                vel = [(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in masses]
                fixed = [False, False, True, True]
                alpha = sum((lx*lx+ly*ly)/m for (lx,ly),m,f in zip(modes,masses,fixed) if not f)
                r0 = sum(lx*vx+ly*vy for (lx,ly),(vx,vy) in zip(modes,vel))
                rf = sum(lx*vx+ly*vy for (lx,ly),(vx,vy),f in zip(modes,vel,fixed) if f)
                s = rng.uniform(1.e-5, 2)
                impulse = r0*(-math.expm1(-s))/alpha
                heat = r0*r0*(-math.expm1(-2*s))/(2*alpha)
                new = [(vx,vy) if f else (vx-impulse*lx/m,vy-impulse*ly/m)
                       for (vx,vy),(lx,ly),m,f in zip(vel,modes,masses,fixed)]
                dK = sum(.5*m*((nx-vx)*(nx+vx)+(ny-vy)*(ny+vy))
                         for m,(vx,vy),(nx,ny) in zip(masses,vel,new))
                work = impulse*rf
                sink = heat if fixed_owner else 0
                dU = 0 if fixed_owner else heat
                with self.subTest(fixed_owner=fixed_owner,case=case):
                    self.assertAlmostEqual(dK+dU,work-sink,delta=2.e-12*(1+abs(work)+heat))
                    for f,old,n in zip(fixed,vel,new):
                        if f: self.assertEqual(old,n)
                    for k in (0,1):
                        dp = sum(m*(nv[k]-ov[k]) for m,ov,nv in zip(masses,vel,new))
                        boundary_dp = impulse*sum(l[k] for l,f in zip(modes,fixed) if f)
                        self.assertAlmostEqual(dp,boundary_dp,delta=2.e-12*(1+abs(boundary_dp)))
    def test_static_boundary_work_is_zero(self):
        alpha,r0,rf,s=2.,3.,0.,.4
        impulse=r0*(-math.expm1(-s))/alpha
        heat=r0*r0*(-math.expm1(-2*s))/(2*alpha)
        dK=-impulse*(r0-rf)+.5*alpha*impulse*impulse
        self.assertAlmostEqual(dK+heat,0.)

if __name__ == '__main__':
    unittest.main()
