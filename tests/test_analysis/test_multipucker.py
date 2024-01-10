#!/usr/bin/env python

import unittest
import pytraj as pt
from pytraj.utils import fn
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

def test_multipucker_furanoid():
    print('MultiPucker 5-member ring pucker, Cremer & Pople Furanoid test')

    cm = '''
    multipucker Furanoid puckertype furanoid:C2:C3:C4:C5:O2 cremer out furanoid.dat amplitude ampout furanoid.dat range360
    '''
    traj = pt.iterload(fn('Test_Pucker/Furanoid.mol2'))
    state = pt.load_cpptraj_state(cm,traj)
    state.run()
    
    data_fura = pt.multipucker(traj=traj, name='Furanoid',
                                puckertype=['furanoid','C2','C3','C4','C5','O2'],
                                method='cremer', out='furanoid.dat',
                                amplitude=True, amp_out='furanoid.dat',
                                range360=True)
    
    aa_eq(data_fura.values, state.data[1:].values)
def test_multipucker_pyranoid():
    print('MultiPucker 6-member ring pucker, Cremer & Pople Pyranoid test')
    
    cm = '''
    multipucker Pyranoid puckertype pyranoid:C1:C2:C3:C4:C5:O5 cremer out pyranoid.type.dat amplitude ampout pyranoid.type.dat theta thetaout pyranoid.type.dat range360
    multipucker Pyr pyranose cremer out pyranoid.auto.dat amplitude ampout pyranoid.auto.dat theta thetaout pyranoid.auto.dat range360
    '''
    traj = pt.iterload(fn('Test_Pucker/Pyranoid.mol2'))
    state = pt.load_cpptraj_state(cm,traj)
    state.run()

    data_pyranoid = pt.multipucker(traj=traj, name='Pyranoid',
                                    puckertype=['pyranoid','C1','C2','C3','C4','C5','05'],
                                    method='cremer', out='pyranoid.type.dat', 
                                    amplitude=True, amp_out='pyranoid.type.dat',
                                    theta= True, theta_out='pyranoid.type.dat',
                                    range360=True)
    data_pyranose = pt.multipucker(traj=traj, name='Pyr',
                                    puckertype='pyranose',
                                    method='cremer', out='pyranoid.auto.dat',
                                    amplitude=True, amp_out='pyranoid.auto.dat',
                                    theta=True, theta_out='pyranoid.auto.dat',
                                    range360=True)
    
    aa_eq(data_pyranoid.values, state.data[1:4])
    aa_eq(data_pyranose.values, state.data[4:])