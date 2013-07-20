import sys
import itertools
from simtk import unit
import simtk.openmm as mm
from simtk.openmm import app
import mdtraj as md
from openmmtools.simulation import AsyncSimulation
from mdtraj.reporters import HDF5Reporter
from ipcfg.progressreporter import ProgressReporter

#----------------------------------
TRAJ = '../250.unfolded.confs.h5'
TOP = '../250.unfolded.confs.h5.0.pdb'
TEMPERATURE = 360*unit.kelvin
TIMESTEP = 2*unit.femtoseconds
N_STEPS = int(25*unit.picoseconds / TIMESTEP)
REPORT_INTERVAL = int(5*unit.picoseconds / TIMESTEP)
N_DEVICES = 1

#----------------------------------
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
    integrator = mm.LangevinIntegrator(TEMPERATURE, 1.0/unit.picoseconds, TIMESTEP)
    integrator.setConstraintTolerance(0.00001)
    simulations.append(AsyncSimulation(topology, system, integrator, mm.Platform.getPlatformByName("CUDA"),
                                       {'CudaPrecision': 'mixed', 'CudaDeviceIndex': str(i)}))


for frame, simulation in zip(range(3), itertools.cycle(simulations)):
    print 'Device Busy? %s' % simulation.isBusy()
    if simulation.isBusy():
        print 'Frame %d, waiting for device' % frame
        simulation.wait()
        print 'Frame %d, done waiting.' % frame
    print 'CUDA Device Index: %s' % simulation.context.getPlatform().getPropertyValue(simulation.context, 'CudaDeviceIndex')
    print "Equilibrating frame %d" % frame

    outfn = 'unfolded.equilibration.NPT.{}.h5'.format(frame)
    positions = t.xyz[frame]
    box = t.unitcell_lengths[frame]
    topology.setUnitCellDimensions(unit.Quantity(box.tolist(), unit.nanometers))
    
    print 'Setting positions...'
    simulation.context.setPositions(positions)
    print 'Setting Velocities...'
    simulation.context.setVelocitiesToTemperature(TEMPERATURE)
    print 'Setting unitcell vectors...'
    simulation.context.setPeriodicBoxVectors([box[0],0,0], [0,box[1],0], [0,0,box[2]])
    print '(Re)setting clock...'
    simulation.context.setTime(0.0)

    # Remove any existing reporters from the previous round
    print 'Adding reporters'
    simulation.reporters.append(HDF5Reporter(outfn, REPORT_INTERVAL))
    simulation.reporters.append(ProgressReporter(sys.stdout, REPORT_INTERVAL, N_STEPS))
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    for r in simulation.reporters:
        r.report(simulation, state)

    def callback():
        print 'Removing reporters for frame %d' % frame
        del simulation.reporters[:]
        
    simulation.asyncstep(N_STEPS, callback=callback)
