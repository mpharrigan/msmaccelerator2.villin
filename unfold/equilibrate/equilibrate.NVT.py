import sys
import numpy as np
import threading
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

TOP = '../250.unfolded.confs.h5.0.pdb'
TEMPERATURE = 360*unit.kelvin
TIMESTEP = 2*unit.femtoseconds
N_STEPS = int(50*unit.picoseconds / TIMESTEP)
REPORT_INTERVAL = int(10*unit.picoseconds / TIMESTEP)
N_FRAMES = 250

# Average box lengths from the NPT equilibrations
BOX = mm.Vec3(1, 1, 1) * 5.4511113281249997 * unit.nanometers

#--------------------------------------------------------

PLATFORM = mm.Platform.getPlatformByName('OpenCL')
#N_DEVICES = int(subprocess.check_output('deviceCount'))
N_DEVICES = 1
DEVICE_INDEX = 'OpenCLDeviceIndex'
PRECISION = 'OpenCLPrecision'

#N_STEPS_PER_DEVICE = int(N_STEPS * N_FRAMES / float(N_DEVICES))

#---------------------------------------------------------
print 'N_DEVICES', N_DEVICES

app.Topology.loadBondDefinitions('../../residues-nle.xml')
topology = app.PDBFile(TOP).topology
topology2 = md.load(TOP).topology

forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml', '../../amber99sbildn-nle.xml')
system = forcefield.createSystem(topology, nonbondedMethod=app.PME, 
                                 nonbondedCutoff=9.5*unit.angstroms,
                                 constraints=app.HBonds, rigidWater=True, 
                                 ewaldErrorTolerance=0.0005)

simulations = []
for i in range(N_DEVICES):
    print 'Initializing device %d' % i
    integrator = mm.LangevinIntegrator(TEMPERATURE, 1.0/unit.picoseconds, TIMESTEP)
    integrator.setConstraintTolerance(0.00001)
    s = AsyncSimulation(topology, system, integrator, PLATFORM, {
        PRECISION: 'Mixed', DEVICE_INDEX: str(i)})
    #s.reporters.append(ProgressReporter("NVT.%d.log" % i, REPORT_INTERVAL, N_STEPS))
    s.reporters.append(ProgressReporter(sys.stdout, REPORT_INTERVAL, N_STEPS))
    simulations.append(s)
    

NVTOUT = md.HDF5TrajectoryFile('NVT.h5', 'w')
NVTOUT.topology = topology
NVTOUTLOCK = threading.Lock()

for frame, simulation in zip(range(N_FRAMES), itertools.cycle(simulations)):
    nptfn = 'unfolded.equilibration.NPT.{}.dcd'.format(frame)
    if simulation.isBusy():
        print 'Frame [%d], starting to wait for device' % frame
        simulation.wait()
        print 'Frame [%d], done waiting for device.' % frame
    print 'Frame [%d], acquired device.' % frame
    print 'Device Index: %s' % PLATFORM.getPropertyValue(simulation.context, DEVICE_INDEX)
    print "Equilibrating [%d]" % frame

    positions = md.load(nptfn, top=topology2).xyz[-1]
    
    print 'Setting positions [%d]...' % frame
    simulation.context.setPositions(positions)
    print 'Setting Velocities [%d]...' % frame
    simulation.context.setVelocitiesToTemperature(TEMPERATURE)
    print 'Setting unitcell vectors [%d]...' % frame
    simulation.context.setPeriodicBoxVectors([BOX[0],0,0], [0,BOX[1],0], [0,0,BOX[2]])
    print '(Re)setting clock [%d]...' % frame
    simulation.context.setTime(0.0)
    simulation.currentStep = 0
    print 'Minimizing a few iterations [%d]...' % frame
    simulation.minimizeEnergy(maxIterations=100)
    
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    for r in simulation.reporters:
        r.report(simulation, state)

    #def write():
    #    NVTOUTLOCK.acquire()
    #    state = simulation.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=True)
    #    coordinates = state.getPositions(asNumpy=True)
    #    vecs = state.getPeriodicBoxVectors(asNumpy=True)
    #    NVTOUT.write(coordinates=coordinates, time=frame,
    #                 cell_lengths=np.diag(vecs),
    #                 cell_angles=np.array([90.0, 90.0, 90.0]))
    #    NVTOUT.flush()
    #    NVTOUTLOCK.release()
                     
    simulation.step(N_STEPS)
    state = simulation.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=True)
    coordinates = state.getPositions(asNumpy=True)
    vecs = state.getPeriodicBoxVectors(asNumpy=True)
    NVTOUT.write(coordinates=coordinates, time=frame,
                 cell_lengths=np.diag(vecs),
                 cell_angles=np.array([90.0, 90.0, 90.0]))


