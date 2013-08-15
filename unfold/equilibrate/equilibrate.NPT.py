import sys
import subprocess
import itertools
from simtk import unit
import simtk.openmm as mm
from simtk.openmm import app
import mdtraj as md
from openmmtools.simulation import AsyncSimulation
from mdtraj.reporters import HDF5Reporter
from ipcfg.progressreporter import ProgressReporter

#--------------------------------------------------------

TRAJ = '../250.unfolded.confs.h5'
TOP = '../250.unfolded.confs.h5.0.pdb'
TEMPERATURE = 360*unit.kelvin
TIMESTEP = 2*unit.femtoseconds
N_STEPS = int(50*unit.picoseconds / TIMESTEP)
REPORT_INTERVAL = int(10*unit.picoseconds / TIMESTEP)
N_FRAMES = 250

#--------------------------------------------------------

PLATFORM = mm.Platform.getPlatformByName('OpenCL')
N_DEVICES = int(subprocess.check_output('deviceCount'))
DEVICE_INDEX = 'OpenCLDeviceIndex'
PRECISION = 'OpenCLPrecision'

#---------------------------------------------------------
print 'N_DEVICES', N_DEVICES

app.Topology.loadBondDefinitions('../../residues-nle.xml')
t = md.load(TRAJ)
topology = app.PDBFile(TOP).topology

forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml', '../../amber99sbildn-nle.xml')
system = forcefield.createSystem(topology, nonbondedMethod=app.PME, 
                                 nonbondedCutoff=9.5*unit.angstroms,
                                 constraints=app.HBonds, rigidWater=True, 
                                 ewaldErrorTolerance=0.0005)
system.addForce(mm.MonteCarloBarostat(1*unit.atmosphere, TEMPERATURE, 25))


simulations = []
for i in range(N_DEVICES):
    print 'Initializing device %d' % i
    integrator = mm.LangevinIntegrator(TEMPERATURE, 1.0/unit.picoseconds, TIMESTEP)
    integrator.setConstraintTolerance(0.00001)
    s = AsyncSimulation(topology, system, integrator, PLATFORM, {
        PRECISION: 'Mixed', DEVICE_INDEX: str(i)})
    simulations.append(s)
    

for frame, simulation in zip(range(N_FRAMES), itertools.cycle(simulations)):
    if simulation.isBusy():
        print 'Frame [%d], starting to wait for device' % frame
        simulation.wait()
        print 'Frame [%d], done waiting for device.' % frame
    print 'Frame [%d], acquired device.' % frame
    print 'Device Index: %s' % PLATFORM.getPropertyValue(simulation.context, DEVICE_INDEX)
    print "Equilibrating [%d]" % frame

    dcdfn = 'unfolded.equilibration.NPT.{}.dcd'.format(frame)
    progressfn = 'unfolded.equilibration.NPT.{}.log'.format(frame)
    positions = t.xyz[frame]
    box = t.unitcell_lengths[frame]
    
    print 'Setting positions [%d]...' % frame
    simulation.context.setPositions(positions)
    print 'Setting Velocities [%d]...' % frame
    simulation.context.setVelocitiesToTemperature(TEMPERATURE)
    print 'Setting unitcell vectors [%d]...' % frame
    simulation.context.setPeriodicBoxVectors([box[0],0,0], [0,box[1],0], [0,0,box[2]])
    print '(Re)setting clock [%d]...' % frame
    simulation.context.setTime(0.0)

    # Remove any existing reporters from the previous round
    del simulation.reporters[:]
    simulation.reporters.append(app.DCDReporter(dcdfn, REPORT_INTERVAL))
    simulation.reporters.append(ProgressReporter(progressfn, REPORT_INTERVAL, N_STEPS))
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    for r in simulation.reporters:
        r.report(simulation, state)

    simulation.asyncstep(N_STEPS)

