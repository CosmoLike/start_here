#!/usr/bin/python
#from __future__ import print_function
# import sys
# sys.path.insert(0, '/usr/local/lib/python2.7/site-packages/')
import matplotlib
from getdist import plots, MCSamples
import getdist, IPython
print('Version: ',getdist.__version__)
import matplotlib.pyplot as plt
print matplotlib.__version__
#Get some random samples for demonstration:
#make random covariance, then independent samples from Gaussian
import numpy as np

zp = 0.2
names = ["Om","s8","ns","w0","wa","Ob","h0"]
labels =["\Omega_m","\sigma_8","n_s","w_0","w_\mathrm{a}","\Omega_b","h_0"]

def single_plot(filename, legends, limits, out):
	chain1 = np.genfromtxt("like/"+filename[0])[:,:7]
	chain2 = np.genfromtxt("like/"+filename[1])[:,:7]
	chain1[:,3] = chain1[:,3]+(1-1./(1+0.266))*chain1[:,4]
	chain2[:,3] = chain2[:,3]+(1-1./(1+zp))*chain2[:,4]
	
	labels[3] ="w_\mathrm{p}"
	names[3] = "wp"

	samples1 = MCSamples(samples=chain1,names = names, labels = labels)
	samples2 = MCSamples(samples=chain2,names = names, labels = labels)
	samples1.fine_bins_2D = 20
	samples2.fine_bins_2D = 20
	
	g = plots.getSinglePlotter()
	g.plot_2d([samples1, samples2],"wp","wa", filled=[False,True],lims=limits,colors=['red','blue'])
	g.settings.legend_fontsize = 12
	g.add_legend(legends)
	g.export('plots/'+out)

limits=[-1.1, -0.9, -0.5, 0.5]
filename=["like_prior_Planck_BAO_SN","like_LSST_Ntomo10_2pt_clusterN_clusterWL_no_sys"]
legends =["Current BOSS+JLA+PL","LSST multi-probe"]
single_plot(filename,legends,limits,"LSST_vs_today.pdf")

