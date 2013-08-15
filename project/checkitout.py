from msmaccelerator.core import markovstatemodel
from msmbuilder import msm_analysis as msma

def l(dir):
	msm = markovstatemodel.MarkovStateModel.load(dir)
	return msm
