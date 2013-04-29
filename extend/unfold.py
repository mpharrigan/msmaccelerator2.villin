import sys
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

import numpy as np
import IPython as ip

##############################################################################

timestep = 2.0*femtoseconds
elongation_factor = 8.0
total_time = 10*nanoseconds
n_steps = int(total_time / timestep)
report_interval = int(10*picoseconds / timestep)
temperature = 700*kelvin

##############################################################################

Topology.loadBondDefinitions('../residues-nle.xml')
pdb = PDBFile('../native.pdb')

forcefield = ForceField('amber99sbildn.xml', '../amber99sbildn-nle.xml',
                        '../amber99-nle_obc.xml')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic, 
    constraints=HBonds, rigidWater=True, nonbondedCutoff=1.0*nanometers)

integrator = LangevinIntegrator(temperature, 91.0/picoseconds, timestep)
integrator.setConstraintTolerance(1e-5)

simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)


print 'Minimizing...'
simulation.minimizeEnergy()

print 'Setting velocities...'
simulation.context.setVelocitiesToTemperature(temperature)

# setup the reporters
simulation.reporters.append(DCDReporter('unfold.py.dcd', report_interval))
simulation.reporters.append(StateDataReporter(sys.stdout, report_interval, temperature=True,
                                              step=True, time=True,
                                              kineticEnergy=True, potentialEnergy=True))    
print 'Starting dynamics!'
simulation.step(n_steps)

positions = simulation.context.getState(getPositions=True).getPositions()
with open('unfold.py.pdb', 'w') as f:
    PDBFile.writeFile(pdb.topology, positions, f)
