"""
Solvate every frame in a molecular dynamics trajectory with the same number of
water molecules in a cubic box.
"""
#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from __future__ import division
import os
import sys
import time

import numpy as np
import scipy.optimize

import mdtraj as md        # Available at https://pypi.python.org/pypi/mdtraj
from simtk import unit
import simtk.openmm as mm
from simtk.openmm import app

#-----------------------------------------------------------------------------
# Globals. These settings can all be altered.
#-----------------------------------------------------------------------------

# Load up custom residue definitions. This is necessary for OpenMM to recognize
# the unnnatural amino acid that I'm using.
app.Topology.loadBondDefinitions('../residues-nle.xml')

PADDING = 1.0*unit.nanometers     # Minimum distance from protein to box edge.
TOPOLOGY = '../native.pdb'        # File with the system topology.
TRAJECTORY = '../extend/unfold.py.dcd'  # Trajectory of frames to be solvated.
STRIDE = 1                        # Load up only every stride-th frame from TRAJECTORY

# ForceField to use for determining van der Waals radii and atomic charges
# when adding solvent.
FF = app.ForceField('amber99sbildn.xml', 'tip3p.xml', '../amber99sbildn-nle.xml')
SITES_PER_WATER = 3

OUT_TOPOLOGY = 'solvated.pdb'             # Solvated topology will be saved here.
SOLVATED_DCD = 'solvated.dcd'  # The solvated systems, before running equilibration.

assert not os.path.exists(OUT_TOPOLOGY), "%s already exists." % OUT_TOPOLOGY
assert not os.path.exists(SOLVATED_DCD), "%s already exists." % SOLVATED_DCD

#-----------------------------------------------------------------------------
# MD settings for minimization / equilibration
#-----------------------------------------------------------------------------

# N_STEPS_MINIMIZATION = 100

# TIMESTEP = 2*unit.femtoseconds
# TEMPERATURE = 300*unit.kelvin
# N_STEPS_EQUILIBRATION = int(5 * unit.picoseconds / TIMESTEP)

# def createSystem(topology):
#     system = FF.createSystem(topology, nonbonedMethod=app.PME,
#                              nonbondedCutoff=9*unit.angstroms,
#                              constraints=app.HBonds)
#     system.addForce(mm.MonteCarloBarostat(1*unit.atmosphere, TEMPERATURE))
#     return system

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

def main():
    traj = md.load(TRAJECTORY, stride=STRIDE, top=TOPOLOGY)
    
    print 'Rotating frames so that they fit in the smallest possible box...'
    for i in progress(range(traj.n_frames)):
        traj.xyz[i] = rotate_to_pack(traj.xyz[i])

    top = app.PDBFile(TOPOLOGY).topology

    #boxes = np.empty((traj.n_frames, 3))
    n_atoms = np.empty(traj.n_frames)

    #print 'Solvating all of the frames with padding %s...' % PADDING
    #for i in progress(range(traj.n_frames)):
    #    modeller = app.Modeller(top, traj.xyz[i].tolist())
    #    modeller.addSolvent(FF, model='tip3p', boxSize=None, padding=PADDING)

    #    boxes[i] = modeller.topology.getUnitCellDimensions().value_in_unit(unit.nanometers)
    #    n_atoms[i] = md.utils.ilen(modeller.topology.atoms())

    # Now, find the largest box, and solvate all of the frames to that box size.
    #print 'Initial number of atoms in each solvated box:'
    #print n_atoms
    #print 'Initial box volumes:'
    #print np.prod(boxes, axis=1)

    #largest_frame = np.argmax(np.prod(boxes, axis=1))
    #largest_box = unit.Quantity(boxes[largest_frame].tolist(), unit.nanometers)
    largest_box = unit.Quantity([8.067393779754639, 8.067393779754639, 8.067393779754639], unit.nanometers)
    #print 'Largest box: %d' % largest_frame
    #print 'Dimensions:  %s' % largest_box
    #print 'N atoms   :  %d' % n_atoms[largest_frame]

    print '\nRe-solvating all frames into this box...'
    modellers = []
    for i in progress(range(traj.n_frames)):
        modeller = app.Modeller(top, traj.xyz[i].tolist())
        modeller.addSolvent(FF, model='tip3p', boxSize=largest_box)
        n_atoms[i] = md.utils.ilen(modeller.topology.atoms())
        modellers.append(modeller)

    print 'After resolvation, the number of atoms in each solvated box:'
    print n_atoms

    n_remove_by_frame = np.asarray((n_atoms - np.min(n_atoms)) / SITES_PER_WATER, dtype=int)
    print '\nNumber of waters to remove from each frame:'
    print n_remove_by_frame
    print '\nRemoving waters from each frame until we get the same number per box'

    for modeller, n_remove in progress(zip(modellers, n_remove_by_frame)):
        waters = [res for res in modeller.topology.residues() if res.name == 'HOH']
        np.random.shuffle(waters)
        modeller.delete(waters[0:n_remove])

    print 'Final number of atoms in each frame'
    print [md.utils.ilen(mod.topology.atoms()) for mod in modellers]

    with open(OUT_TOPOLOGY, 'w') as f:
        # Debugging an issue where there seems to be a confusion between
        # positions having some floats in it and some np.float64s. Lets try to
        # standardize it.
        positions = unit.Quantity(np.array(modellers[0].positions._value).tolist(), unit.nanometers)
        app.PDBFile.writeFile(modellers[0].topology, positions, f)

    with open(SOLVATED_DCD, 'w') as f:
        print 'Saving solvated structures'
        dcdfile = app.DCDFile(f, modellers[0].topology, TIMESTEP)
        for modeller in progress(modellers):
            assert all(a.name == b.name for a, b in zip(modeller.topology.atoms(), modellers[0].topology.atoms()))
            dcdfile.writeModel(modeller.positions, modeller.topology.getUnitCellDimensions())


def rotate_to_pack(points):
    """Rotate a set of points to fit in the smallest possible axis-aligned
    box.

    For example, if the points are roughly linear, we'd rotate them to
    point along the (1, 1, 1) vector so that they'd fall along the diagonal
    of an axis-aligned box.

    Parameters
    ----------
    points : np.ndarray, shape=(N, 3)
        A set of points in 3D.

    Returns
    -------
    result : np.ndarray, shape=(N, 3)
        A rotated a translated version of `points`.
    """
    points = np.asarray(points, dtype=float)
    assert points.ndim == 2 and points.shape[1] == 3, 'Shape Error'
    points = points - np.mean(points, axis=0)

    def rotation_matrix(angles):
        """Calculate a 3x3 rotation matrix from Euler angles using
        the sxyz axis sequence.

        This code is adapted from transformations.py, copyright
        2006-2013, Christoph Gohlke, and released under the 3-clause BSD.
        http://www.lfd.uci.edu/~gohlke/code/transformations.py.html.
        """
        assert len(angles) == 3
        si, sj, sk = np.sin(angles)
        ci, cj, ck = np.cos(angles)
        cc, cs = ci*ck, ci*sk
        sc, ss = si*ck, si*sk

        M = np.array([[cj*ck, sj*sc-cs, sj*cc+ss],
                      [cj*sk, sj*ss+cc, sj*cs-sc],
                      [-sj, cj*si, cj*ci]])
        return M

    def objective(euler_angles):
        """The maximum axis-aligned width of the points, under a specified
        rotation"""

        result = np.dot(points, rotation_matrix(euler_angles))
        width = np.max(np.max(result, axis=0) - np.min(result, axis=0))
        return width

    euler_angles = scipy.optimize.fmin(objective, [0.0, 0.0, 0.0], disp=False)
    return np.dot(points, rotation_matrix(euler_angles))


def progress(iterable, max=None):
    "Wrapper around an iterable that prints its progress as its consumed"
    i = 0
    if hasattr(iterable, '__len__'):
        max = len(iterable)
    if max is None:
        max = '?'

    start = time.time()
    for item in iterable:
        speed = (time.time() - start) / i if i > 0 else np.nan
        if max:
            sys.stdout.write('\r%d / %s (%.5f s/iter) ' % (i+1, max, speed))

        sys.stdout.flush()
        i += 1

        yield item

    print '\n'

if __name__ == '__main__':
    main()
