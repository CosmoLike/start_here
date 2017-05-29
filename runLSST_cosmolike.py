#!/usr/bin/env /Users/teifler/Python-2.7.8/python.exe 

import os

from cosmolike_libs import * 

file_source_z = os.path.join(dirname, "../../zdistris/zdistribution_LSST")
file_lens_z = os.path.join(dirname, "../../zdistris/zdistribution_const_comoving")
data_file = os.path.join(dirname, "datav/LSST_all_2pt_clusterN_clusterWL_fid")
cov_file = os.path.join(dirname, "cov/cov_LSST_2.600000e+01_1.800000e+04_Rmin10_Ncl25_Ntomo10_2pt_clusterN_clusterWL_inv")
chain_file = os.path.join(dirname, "./like/like_LSST_Ntomo10_2pt_clusterN_clusterWL_no_sys")

initcosmo()
initbins(25,20.0,5000.0,5000.0,10.0,7)
initsurvey("LSST")
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","redmagic")
initclusters()
initia("none","DEEP2")
initpriors("none","none","PhotoBAO","none")
initprobes("all_2pt_clusterN_clusterWL")
initdatainv(cov_file,data_file)

sample_params= sample_cosmology_only()
#sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear())
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 
sample_main(sample_params,4000,14,8,chain_file, blind=False)
