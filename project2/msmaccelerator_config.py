# Configuration file for msmaccelerator.

c = get_config()
from simtk.unit import angstrom, nanometer, picoseconds, femtoseconds, nanoseconds

#------------------------------------------------------------------------------
# App configuration
#------------------------------------------------------------------------------

# This is an application.

# App will inherit config from: Application

# The Logging format template
# c.App.log_format = '[%(name)s] %(message)s'

# Set the log level by value or name.
# c.App.log_level = 30

#------------------------------------------------------------------------------
# RootApplication configuration
#------------------------------------------------------------------------------

# This is an application.

# RootApplication will inherit config from: App, Application

# The Logging format template
# c.RootApplication.log_format = '[%(name)s] %(message)s'

# Set the log level by value or name.
# c.RootApplication.log_level = 30

#------------------------------------------------------------------------------
# Device configuration
#------------------------------------------------------------------------------

# This is an application.

# Device will inherit config from: App, Application

# ZeroMQ port to connect to the server on
#c.Device.zmq_port = 12347

# The Logging format template
# c.Device.log_format = '[%(name)s] %(message)s'

# Set the log level by value or name.
# c.Device.log_level = 30

# URL to connect to server with
c.Device.zmq_url = 'icmemaster'

#------------------------------------------------------------------------------
# Modeler configuration
#------------------------------------------------------------------------------

# This is an application.

# Modeler will inherit config from: Device, App, Application

# URL to connect to server with
# c.Modeler.zmq_url = '127.0.0.1'

# Set the log level by value or name.
# c.Modeler.log_level = 30

# Should we use a custom distance metric for clusering instead of RMSD?
# c.Modeler.use_custom_metric = False

# The method used for clustering structures in the MSM.
c.Modeler.clusterer = 'hybrid'

# Symmetrization method for constructing the reversibile counts matrix.
# c.Modeler.symmetrize = None

# File containing a pickled metric for use in clustering.
# c.Modeler.custom_metric_path = 'metric.pickl'

# Subsample data by taking only every stride-th point
# c.Modeler.stride = 1

# Lag time for building the model, in units of the stride. Currently, we are not
# doing the step in MSMBuilder that is refered to as "assignment", where you
# assign the remaining data that was not used during clustering to the cluster
# centers that were identified.
c.Modeler.lag_time = 100

# File containing the indices of atoms to use in the RMSD computation. Using a
# PDB as input, this file can be created with the MSMBuilder script
# CreateAtomIndices.py

c.Modeler.rmsd_atom_indices = 'AtomIndices.dat'

# Do ergodic trimming when constructing the Markov state model. This is
# generally a good idea for building MSMs in the high-data regime where you wish
# to prevent transitions that appear nonergodic because they've been
# undersampled from influencing your model, but is inappropriate in the sparse-
# data regime when you're using min-counts sampling, because these are
# precisiely the states that you're most interested in.
# c.Modeler.ergodic_trimming = False

# PDB file giving the topology of the system
# c.Modeler.topology_pdb = <Undefined>

# ZeroMQ port to connect to the server on
# c.Modeler.zmq_port = 12345

# The Logging format template
# c.Modeler.log_format = '[%(name)s] %(message)s'

# Distance cutoff for clustering, in nanometers. We will continue to create new
# clusters until each data point is within this cutoff from its cluster center.

c.Modeler.clustering_distance_cutoff = 0.5

#------------------------------------------------------------------------------
# OpenMMSimulator configuration
#------------------------------------------------------------------------------

# This is an application.

# OpenMMSimulator will inherit config from: Device, App, Application

# Path to the XML file containing the OpenMM Integrator to use
# c.OpenMMSimulator.integrator_xml = 'integrator.xml'

# Set the log level by value or name.
# c.OpenMMSimulator.log_level = 30

# Interval at which to save positions to a disk, in units of steps
# c.OpenMMSimulator.report_interval = 1000

# Path to the XML file containing the OpenMM system to propagate
# c.OpenMMSimulator.system_xml = 'system.xml'

# Do local energy minimization on the configuration that's passed to me, before
# running dynamics
# c.OpenMMSimulator.minimize = True

# OpenMM device index for CUDA or OpenCL platforms. This is used to select which
# GPU will be used on a multi-gpu system. This option is ignored on reference
# platform
# c.OpenMMSimulator.device_index = 0

# URL to connect to server with
# c.OpenMMSimulator.zmq_url = '127.0.0.1'

# The OpenMM platform on which to run the simulation
# c.OpenMMSimulator.platform = 'CUDA'

# Number of steps of dynamics to do
# c.OpenMMSimulator.number_of_steps = 10000

# ZeroMQ port to connect to the server on
# c.OpenMMSimulator.zmq_port = 12345

# Choose random initial velocities from the Maxwell-Boltzmann distribution
# c.OpenMMSimulator.random_initial_velocities = True

# The Logging format template
# c.OpenMMSimulator.log_format = '[%(name)s] %(message)s'

# Number of steps of dynamics to do
c.OpenMMSimulator.number_of_steps = int((15 * nanoseconds) / (2.0 * femtoseconds))

# Interval at which to save positions to a disk, in units of steps
c.OpenMMSimulator.report_interval = int((100 * picoseconds) / (2.0 * femtoseconds))

#------------------------------------------------------------------------------
# AmberSimulator configuration
#------------------------------------------------------------------------------

# This is an application.

# AmberSimulator will inherit config from: Device, App, Application

# Which AMBER executable to use?
# c.AmberSimulator.executable = 'pmemd'

# Directory to work in. If not set, we'll requirest a temporary directory from
# the OS and clean it up when we're finished. This option is useful for
# debugging.
# c.AmberSimulator.workdir = ''

# Set the log level by value or name.
# c.AmberSimulator.log_level = 30

# Parameter/topology file for the system
# c.AmberSimulator.prmtop = <Undefined>

# AMBER .in file controlling the production run. If no production is desired, do
# not set this parameter.
# c.AmberSimulator.mdin = <Undefined>

# URL to connect to server with
# c.AmberSimulator.zmq_url = '127.0.0.1'

# The Logging format template
# c.AmberSimulator.log_format = '[%(name)s] %(message)s'

# ZeroMQ port to connect to the server on
# c.AmberSimulator.zmq_port = 12345

#------------------------------------------------------------------------------
# Interactor configuration
#------------------------------------------------------------------------------

# This is an application.

# Interactor will inherit config from: Device, App, Application

# Go into interactive shell mode, in which you can interact with the server via
# an IPython read-eval-print loop.
# c.Interactor.shell = False

# Set the log level by value or name.
# c.Interactor.log_level = 30

# URL to connect to server with
# c.Interactor.zmq_url = '127.0.0.1'

# Set the server's beta parameter on the fly, changing the balance between
# exploration and optimization of the already discovered pathways/rates. For
# details this parameter affects the sampling. See the help text for the
# CountsSampler.beta trait
# c.Interactor.set_beta = None

# The Logging format template
# c.Interactor.log_format = '[%(name)s] %(message)s'

# ZeroMQ port to connect to the server on
# c.Interactor.zmq_port = 12345

# Path to the database (sqlite3 file)
# c.Interactor.db_path = 'db.sqlite'

#------------------------------------------------------------------------------
# BaseServer configuration
#------------------------------------------------------------------------------

# This is an application.

# BaseServer will inherit config from: App, Application

# The Logging format template
# c.BaseServer.log_format = '[%(name)s] %(message)s'

# Set the log level by value or name.
# c.BaseServer.log_level = 30

# ZeroMQ port to serve on
# c.BaseServer.zmq_port = 12345

# Path to the database (sqlite3 file)
# c.BaseServer.db_path = 'db.sqlite'

#------------------------------------------------------------------------------
# AdaptiveServer configuration
#------------------------------------------------------------------------------

# This is an application.

# AdaptiveServer will inherit config from: BaseServer, App, Application

# Which MD engine do you want to configure the server to iterface with? If
# 'OpenMM', the server will emit xml-serialized states to simulators that
# connect. If 'AMBER', the server will instead emit inpcrd files to the
# simulators.
# c.AdaptiveServer.md_engine = 'OpenMM'

# Set the log level by value or name.
# c.AdaptiveServer.log_level = 30

# Path to the XML file containing the OpenMM system to propagate. This is
# required by the server, iff md_engine=='OpenMM', to properly serialize the
# starting conformations.
# c.AdaptiveServer.system_xml = 'system.xml'

# Directory on the local filesystem where output trajectories will be saved
# c.AdaptiveServer.traj_outdir = 'trajs/'

# Directory on the local filesystem where MSMs will be saved.
# c.AdaptiveServer.models_outdir = 'models/'

# A PDB used to determine the system's topology. This is sent directly to the
# Simulator. Honestly, I'm not sure exactly why we need it. TODO: ask Peter
# about this.
c.AdaptiveServer.topology_pdb = 'NVT.pdb'

# ZeroMQ port to serve on
# c.AdaptiveServer.zmq_port = 12345

# 'Directory on the local filesystem where starting structures will be stored
# c.AdaptiveServer.starting_states_outdir = 'starting_states'

# The Logging format template
# c.AdaptiveServer.log_format = '[%(name)s] %(message)s'

# Path to the database (sqlite3 file)
# c.AdaptiveServer.db_path = 'db.sqlite'

#------------------------------------------------------------------------------
# MKProfile configuration
#------------------------------------------------------------------------------

# This is an application.

# MKProfile will inherit config from: App, Application

# The Logging format template
# c.MKProfile.log_format = '[%(name)s] %(message)s'

# Set the log level by value or name.
# c.MKProfile.log_level = 30

# Output directory in which to save the file msmbuilder_config.py
# c.MKProfile.output_dir = '.'

#------------------------------------------------------------------------------
# BaseSampler configuration
#------------------------------------------------------------------------------

# Adaptive sampler that simply returns a given initial structure.
#
# This can be used for the very first round when you have no data, and also as a
# base class for the other samplers

# Trajectory file giving the initial structures that you want to sample from.
# This should be a single PDB or other type of loadable trajectory file. These
# structures will only be used in the beginning, before we have an actual MSM to
# use.
c.BaseSampler.seed_structures = 'NVT.h5'

#------------------------------------------------------------------------------
# CentroidSampler configuration
#------------------------------------------------------------------------------

# Adaptive sampler that uses the centroids/"generators"" of the states to choose
# amongst, with a multinomial distibition.
#
# This class DOES NOT actually contain a method to *set* the weights. That is
# done by subclasses. See, for example, CountsSampler, that sets the weights
# using the counts.

# CentroidSampler will inherit config from: BaseSampler

# Trajectory file giving the initial structures that you want to sample from.
# This should be a single PDB or other type of loadable trajectory file. These
# structures will only be used in the beginning, before we have an actual MSM to
# use.
# c.CentroidSampler.seed_structures = 'ala5.pdb'

#------------------------------------------------------------------------------
# CountsSampler configuration
#------------------------------------------------------------------------------

# CountsSampler will inherit config from: CentroidSampler, BaseSampler

# Temperature factor that controls the level of exploration vs. refinement. When
# beta = 0 (high temp), we do full exploration, putting simutions where few
# counts have been seen. When beta = 1 (room temp), we do uniform sampling from
# the microstate, with no preference based on their counts. When beta > 1 (low
# temp), we do refinement, such that we focus on microstates with a high number
# of counts. At beta = 2, we choose microstates proportional to our estimate of
# their current equilibrium propbability. The explicit formula used is: Prob(
# choose state i ) ~ \sum_j C_{ij} ^{ beta - 1 }
# c.CountsSampler.beta = 1

# Trajectory file giving the initial structures that you want to sample from.
# This should be a single PDB or other type of loadable trajectory file. These
# structures will only be used in the beginning, before we have an actual MSM to
# use.
# c.CountsSampler.seed_structures = 'ala5.pdb'

