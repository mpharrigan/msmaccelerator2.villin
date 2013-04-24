import numpy as np
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

###########################
temperature = 300*kelvin
timestep = 1.0*femtosecond
###########################

pdb = PDBFile('solvated.pdb')
pdb.topology.loadBondDefinitions('../residues-nle.xml')
pdb.topology.createStandardBonds()
forcefield = ForceField('amber99sbildn.xml', '../amber99sbildn-nle.xml',
                        'tip3p.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                                 constraints=None, rigidWater=False,
                                 nonbondedCutoff=0.8*nanometers)
integrator = LangevinIntegrator(temperature, 1.0/picosecond, timestep)
integrator.setConstraintTolerance(0.00001)


simulation = Simulation(pdb.topology, system, integrator)# Platform.getPlatformByName('Reference'))
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperature)

state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True,
                                    getParameters=True, enforcePeriodicBox=True)

print 'Kinetic Energy', state.getKineticEnergy()
print 'RMS Velocity', np.sqrt(np.mean((np.square(state.getVelocities()))))

for i in range(10000):
    print 'minimiziation', i
    print simulation.context.getState(getEnergy=True).getPotentialEnergy()
    simulation.minimizeEnergy(maxIterations=50)


positions = simulation.context.getState(getPositions=True).getPositions()
with open('minimized.pdb', 'w') as f:
    PDBFile.writeFile(pdb.topology, positions, f)

#simulation.reporters.append(StateDataReporter(sys.stdout, report_interval, time=True, step=True,
#                                              temperature=True, kineticEnergy=True, potentialEnergy=True))
#simulation.step(2)

#print 'total n_steps. starting.'
#simulation.reporters.append(DCDReporter('equilibration.dcd', report_interval))
#

#print 'Done'
