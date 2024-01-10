#!/usr/bin/env python

import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq

cm = '''
multipucker ADHas resrange 1-3 altona out nucleic.dat
multipucker ADHcp resrange 1-3 cremer out nucleic.dat
'''

def test_multipucker():
    print('MultiPucker command test')
    
    cm = '''
    multipucker ADHas resrange 1-3 altona out nucleic.dat
    multipucker ADHcp resrange 1-3 cremer out nucleic.dat
    '''
    traj = pt.iterload(fn('Test_NAstruct/adh026.3.pdb'))
    state = pt.load_cpptraj_state(cm, traj)
    state.run()
    data_altona = pt.multipucker(traj=traj, out='nucleic.dat', resrange=(1,3),
                                    method='altona')
    data_cremer = pt.multipucker(traj=traj, out='nucleic.dat', resrange=(1,3),
                                    method='cremer')

    aa_eq(data_altona.values, state.data[1:4].values)
    aa_eq(data_cremer.values, state.data[4:].values)