"""Axis extraction uses the same local harmonic convention as the maps."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from plot_2d_fields import FieldMap
from plot_double_core_axis_profiles import positive_axis_profiles


class HarmonicProfiles(unittest.TestCase):
    def test_nonwinding_single_harmonic(self):
        coords = np.array([-2., 0., 2.])
        a = np.zeros((3, 3, 3, 3), complex)
        a[2, 0] = 1 / np.sqrt(2)
        a[2, 1] = -1j / np.sqrt(2)
        field = FieldMap(coords, coords, a, np.zeros((3,3)), np.zeros((3,3,3)), {})
        for profile in positive_axis_profiles(field, 'harmonic'):
            expected = np.zeros((3,3,2), complex)
            expected[1,2] = 1
            np.testing.assert_allclose(profile.values, expected, atol=1e-14)
            np.testing.assert_array_equal(profile.coordinate, [0,2])
        cartesian = positive_axis_profiles(field)
        np.testing.assert_allclose(cartesian[0].values, a[:,:,1,1:])

    def test_bulk_phase_retained(self):
        coords = np.array([-2., 0., 2.])
        a = np.zeros((3,3,3,3), complex)
        x,y = np.meshgrid(coords, coords)
        phase = np.exp(1j*np.arctan2(y,x))
        for i in range(3):
            a[i,i] = phase
        field = FieldMap(coords, coords, a, np.zeros((3,3)), np.zeros((3,3,3)), {})
        xp,yp=positive_axis_profiles(field, 'harmonic')
        self.assertAlmostEqual(xp.values[1,1,-1], 1)
        self.assertAlmostEqual(yp.values[1,1,-1], 1j)


if __name__ == '__main__':
    unittest.main()
