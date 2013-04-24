Villin-NLE-NLE
==============

Files contained herein:

- The gromacs files for amber99sbildn-nle are courtesy of Kyle A. Beauchamp,
  and in the `amber99sb-ildn-nle/` directory.
- `amber99sbildn-nle.xml` contains the forcefield parameters in OpenMM format
  for the NLE residue, translated from the gromacs parameter files.
- `amber99sb-nle_obc.xml` contains the GBSAOBC implicit solvent parameters
  for the NLE residue, translated from the gromacs parameter files.
- `residues-nle.xml` contains the topology description of the NLE residue
  for OpenMM, which is used by the PDBReader to generate the topology.
- `native.pdb` contains the crystal structure, sans water, courtesy of
  Kyle A. Beauchamp.

To get the topology read correctly from the PDB, you need to do something
like this, otherwise the bonds will not be set correctly in the topology.

```
pdb = PDBFile('native.pdb')
pdb.topology.loadBondDefinitions('residues-nle.xml')
pdb.topology.createStandardBonds()
```
  