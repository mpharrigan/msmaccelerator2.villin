from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

###########################
temperature = 300*kelvin
timestep = 2.0*femtosecond
n_steps = int(100*picoseconds / timestep)
report_interval = int(1*picoseconds / timestep)
###########################

Topology.loadBondDefinitions('../residues-nle.xml')
pdb = PDBFile('unfold.py.pdb')

modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField('amber99sbildn.xml', '../amber99sbildn-nle.xml', 'tip3p.xml')

modeller.addHydrogens(forcefield, pH=7.0)
modeller.addSolvent(forcefield, padding=1.0*nanometers)

print 'Creating system...'
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
                                                                  constraints=HBonds, rigidWater=False,
                                                                  nonbondedCutoff=0.8*nanometers)
integrator = LangevinIntegrator(temperature, 1.0/picosecond, timestep)
integrator.setConstraintTolerance(1e-5)

print 'Setting context...'
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
print 'Minimizing...'
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(300*kelvin)
simulation.reporters.append(StateDataReporter(sys.stdout, report_interval, time=True, step=True,
			    temperature=True, kineticEnergy=True, potentialEnergy=True))
simulation.reporters.append(DCDReporter('solvate.py.dcd', report_interval))

print 'Starting dynamics!'
simulation.step(n_steps)

positions = simulation.context.getState(getPositions=True).getPositions()
with open('solvate.py.pdb', 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
        print 'Done'

