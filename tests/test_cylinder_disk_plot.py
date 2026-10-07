import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from plot_cylinder_disk import interpolate,reconstruct,cartesian_to_harmonics

class CylinderDiskTests(unittest.TestCase):
    def test_quadratic_and_nodes(self):
        r=np.array([0.,.2,1.,3.,5.]); q=np.linspace(0,5,101)
        np.testing.assert_allclose(interpolate(r,r*r,q),q*q,atol=1e-13)
        np.testing.assert_allclose(interpolate(r,r*r,r),r*r,atol=1e-13)

    def test_phase_vector_and_mask(self):
        r=np.linspace(0,20,201)
        rng=np.random.default_rng(3)
        a=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        f=reconstruct(r,np.tile(a,(len(r),1,1)),np.tile([1.,2.,3.],(len(r),1)),
                      np.array([3.,21.]),np.array([4.]),1)
        m=np.array([1,0,-1])
        expect=cartesian_to_harmonics(a)*np.exp(1j*(1-m[:,None]-m[None,:])*np.arctan2(4,3))
        np.testing.assert_allclose(cartesian_to_harmonics(f.order_parameter)[:,:,0,0],expect,atol=1e-13)
        np.testing.assert_allclose(f.current[:,0,0],[-1,2,3],atol=1e-13)
        self.assertTrue(np.isnan(f.order_parameter[:,:,0,1]).all())

if __name__=='__main__': unittest.main()
