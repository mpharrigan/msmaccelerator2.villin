##########################################################################
# this script was generated by openmm-builder. to customize it further,
# you can save the file to disk and edit it with your favorite editor.
##########################################################################

#from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

import xml.etree.ElementTree as etree
from openmmtools.pullingforcewrapper import PullingForceWrapper
from openmmtools.webreporter import WebReporter
import threading

import numpy as np
import IPython as ip

##############################################################################

timestep = 2.0*femtoseconds
elongation_factor = 8.0
n_steps = 5000000
n_intervals = 100

##############################################################################

pdb = PDBFile('../native.pdb')

pdb.topology.loadBondDefinitions('../residues-nle.xml')
pdb.topology.createStandardBonds()

forcefield = ForceField('../amber99sbildn.xml', '../amber99sbildn-nle.xml',
                        '../amber99-nle_obc.xml')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic, 
    constraints=None, rigidWater=True, nonbondedCutoff=1.0*nanometers)

integrator = LangevinIntegrator(300*kelvin, 91.0/picoseconds, timestep)
integrator.setConstraintTolerance(0.0001)

pullingforce = PullingForceWrapper(pdb=pdb)
pullingforce.add_to_system(system)

platform = Platform.getPlatformByName('Reference')
simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)


print 'Minimizing...'
simulation.minimizeEnergy()


simulation.context.setVelocitiesToTemperature(300*kelvin)
print 'Equilibrating...'
simulation.step(100)


# setup the reporters
simulation.reporters.append(DCDReporter('pulling.dcd', 10000))
simulation.reporters.append(StateDataReporter(stdout, 10000, temperature=True, step=True, time=True))    


total_time = n_steps * timestep
print 'Going for', total_time.in_units_of(picoseconds)
initial_r0 = pullingforce.get_r0()
final_r0 = elongation_factor * pullingforce.get_r0()

print 'PULLING RATE', ((final_r0 - initial_r0) / total_time).in_units_of(nanometer/nanosecond)


for i, r0 in enumerate(np.linspace(initial_r0, final_r0, n_intervals)):
    print 'moving the pulling points out'
    pullingforce.set_r0(r0, simulation.context)
    simulation.step(n_steps / n_intervals)
