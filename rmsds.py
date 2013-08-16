from __future__ import division
from matplotlib import pyplot as pp
import mdtraj
import numpy as np
import os
import sys



def get_rmsds(traj_fn, native_fn='native.pdb', aindices_fn="project/AtomIndices.dat"):

    aindices = np.loadtxt(aindices_fn, dtype='int')
    native = mdtraj.load(native_fn, atom_indices=aindices)
    traj = mdtraj.load(traj_fn, atom_indices=aindices)

    native = mdtraj.rmsd_cache(native)
    traj = mdtraj.rmsd_cache(traj)

    rmsds = traj.rmsds_to(native, 0)
    return rmsds

def do_all_trajectories(directory='project/trajs/'):

    rr_list = list()
    files = os.listdir(directory)

    for i in xrange(len(files)):
        f = files[i]
        if f.endswith('.h5'):
            # Trajectory

            rr = get_rmsds("%s/%s" % (directory, f))
            rr_list.append(rr)

        print("Did a trajectory! %.2f%% Complete" % (100 * (i / len(files))))


    return rr_list

def plot_some(rr_list, beg=0, number=10, stride=1):
    for rr in rr_list[beg:beg + stride * number:stride]:
        pp.plot(rr)
    pp.ylim(0, 1.5)
    pp.show()

def find_min(rr_list):
    min = sys.maxint
    arg = (-1, -1)
    i = 0
    for rr in rr_list:
        new_min = np.min(rr)
        if new_min < min:
            min = new_min
            arg = (i, np.argmin(rr))
        i += 1

    return min, arg
