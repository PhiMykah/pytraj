from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal, eq_coords
from pytraj import set_world_silent, set_error_silent

fname = "./data/tz2.pdb"


class Test(unittest.TestCase):
    def test_0(self):
        top = Topology("./data/saxs_test/test.pdb")

    def test_1(self):
        top = Topology("./data/saxs_test/test.pdb")

    def test_2_load_pdb(self):
        pdb = mdio.load_pdb(fname)
        traj = mdio.load(fname, fname)
        eq_coords(pdb, traj)


if __name__ == "__main__":
    unittest.main()