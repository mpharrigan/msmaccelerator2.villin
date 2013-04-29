from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

Topology.loadBondDefinitions('../residues-nle.xml')
pdb = PDBFile('solvate.py.pdb')
forcefield = ForceField('amber99sbildn.xml', '../amber99sbildn-nle.xml', 'tip3p.xml')

# create a system
system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic,
                                 nonbondedCutoff=0.8*nanometer, constraints=HBonds,
                                 ewaldErrorTolerance=0.0005)
# and serialize it
print 'writing system.xml file...'
with open('system.xml','w') as f:
    f.write(XmlSerializer.serialize(system))

# do the same for an integrator
integrator = LangevinIntegrator(300*kelvin, 1.0/picoseconds, 2.0*femtoseconds)
integrator.setConstraintTolerance(1e-5)
print 'writing integrator.xml file...'
with open('integrator.xml', 'w') as f:
    f.write(XmlSerializer.serialize(integrator))
