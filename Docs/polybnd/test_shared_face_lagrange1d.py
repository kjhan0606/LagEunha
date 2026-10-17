import unittest
import numpy as np
from shared_face_lagrange1d import initial,step,power_weights,decode

class SharedFaceTests(unittest.TestCase):
    def test_contact_single_step(self):
        y,m=initial(8,'contact',.7);z,ledger,d=step(y,m,.001,(.7,.7))
        np.testing.assert_allclose(z[0],y[0]+.0007,atol=2e-16)
        np.testing.assert_allclose(decode(z,m)[3],.6,atol=5e-15)
        self.assertLess(d['gcl'],2e-16)
    def test_power_envelope(self):
        y,m=initial(8,'acoustic');z,_,_=step(y,m,.002,(0.,0.));w,_=power_weights(z[0],z[1])
        for i in range(8):
            for t in (.1,.5,.9):
                point=(1-t)*z[0][i]+t*z[0][i+1]
                self.assertEqual(int(np.argmin((point-z[1])**2-w)),i)
    def test_energy_ledger(self):
        y,m=initial(8,'compression');z,dl,d=step(y,m,.001,(.1,-.1))
        self.assertLess(abs(sum(z[3])-sum(y[3])-dl[1]),3e-16)
        self.assertGreater(dl[1],0.)
if __name__=='__main__':unittest.main()
