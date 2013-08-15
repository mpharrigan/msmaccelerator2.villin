from __future__ import division
import os
import sys

from simtk import unit
import simtk.openmm as mm
from simtk.openmm import app, XmlSerializer


#-----------------------------------------------------------------------------

NATIVE = '../project/NVT.pdb'
TIMESTEP = 2 * unit.femtoseconds
REPORT_INTERVAL = int((100 * unit.picoseconds) / (2.0 * unit.femtoseconds))
N_STEPS = int((10 * unit.nanoseconds) / (2.0 * unit.femtoseconds))
TEMPERATURE = 360 * unit.kelvin

#-----------------------------------------------------------------------------

# Create System

# Load up custom residue definitions. This is necessary for OpenMM to recognize
# the unnnatural amino acid that I'm using.
app.Topology.loadBondDefinitions('../residues-nle.xml')
pdb = app.PDBFile(NATIVE)
forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml',
                            '../amber99sbildn-nle.xml')



system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME,
                                 nonbondedCutoff=9.5 * unit.angstroms,
                                 constraints=app.HBonds, rigidWater=True,
                                 ewaldErrorTolerance=0.0005)
# system.addForce(mm.MonteCarloBarostat(1 * unit.atmosphere, TEMPERATURE, 25))

print("Writing system.xml")
with open('system.xml', 'w') as f:
    f.write(XmlSerializer.serialize(system))


# Create Integrator

integrator = mm.LangevinIntegrator(TEMPERATURE, 1.0 / unit.picoseconds, TIMESTEP)
integrator.setConstraintTolerance(0.00001)

print("Writing integrator.xml")
with open('integrator.xml', 'w') as f:
    f.write(XmlSerializer.serialize(integrator))

# platform = mm.Platform.getPlatformByName('CUDA')
# properties = {'CudaPrecision': 'mixed'}
# simulation = app.Simulation(modeller.topology, system, integrator, platform,
#                             properties)
# simulation.context.setPositions(modeller.positions)
#
# print 'Setting Velocities'
# simulation.context.setVelocitiesToTemperature(TEMPERATURE)
#
# print 'Adding reporters'
# simulation.reporters.append(HDF5Reporter(OUT_FN, REPORT_INTERVAL))
# simulation.reporters.append(ProgressReporter(sys.stdout, REPORT_INTERVAL, N_STEPS))
#
# state = simulation.context.getState(getPositions=True, getEnergy=True)
# for r in simulation.reporters:
#     r.report(simulation, state)
#
# simulation.step(N_STEPS)
