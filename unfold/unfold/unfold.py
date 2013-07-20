from __future__ import division
import os
import sys

from simtk import unit
import simtk.openmm as mm
from simtk.openmm import app

from mdtraj.reporters import HDF5Reporter
from ipcfg.progressreporter import ProgressReporter

#-----------------------------------------------------------------------------

BOX_SIZE = [54.0, 54.0, 54.0] * unit.angstroms
NATIVE = '../native.pdb'
OUT_FN = 'unfold.tip3p.h5'
TIMESTEP = 2*unit.femtoseconds
REPORT_INTERVAL = int(10*unit.picoseconds / TIMESTEP)
N_STEPS = int(100*unit.nanoseconds / TIMESTEP)
TEMPERATURE = 500*unit.kelvin

assert not os.path.exists(OUT_FN), '%s already exists' % OUT_FN
#-----------------------------------------------------------------------------

# Load up custom residue definitions. This is necessary for OpenMM to recognize
# the unnnatural amino acid that I'm using.
app.Topology.loadBondDefinitions('../residues-nle.xml')
pdb = app.PDBFile(NATIVE)
forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml', '../amber99sbildn-nle.xml')

modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, model='tip3p', boxSize=BOX_SIZE, ionicStrength=40*unit.millimolar)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, 
                                 nonbondedCutoff=9.5*unit.angstroms,
                                 constraints=app.HBonds, rigidWater=True, 
                                 ewaldErrorTolerance=0.0005)
system.addForce(mm.MonteCarloBarostat(1*unit.atmosphere, TEMPERATURE, 25))

integrator = mm.LangevinIntegrator(TEMPERATURE, 1.0/unit.picoseconds, TIMESTEP)
integrator.setConstraintTolerance(0.00001)

platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)

print 'Setting Velocities'
simulation.context.setVelocitiesToTemperature(TEMPERATURE)

print 'Adding reporters'
simulation.reporters.append(HDF5Reporter(OUT_FN, REPORT_INTERVAL))
simulation.reporters.append(ProgressReporter(sys.stdout, REPORT_INTERVAL, N_STEPS))

state = simulation.context.getState(getPositions=True, getEnergy=True)
for r in simulation.reporters:
    r.report(simulation, state)

simulation.step(N_STEPS)
