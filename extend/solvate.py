from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

pdb = PDBFile('unfolded.pdb')

pdb.topology.loadBondDefinitions('../residues-nle.xml')
pdb.topology.createStandardBonds()

modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField('amber99sbildn.xml', '../amber99sbildn-nle.xml', 'tip3p.xml')
#                        '../amber99-nle_obc.xml')

modeller.addHydrogens(forcefield, pH=7.0)
modeller.addSolvent(forcefield, padding=1.0*nanometers)

with open('solvated.pdb', 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print 'Done'
