#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 09:09:03 2019

@author: Sijie
"""

import sys
wd = ('/Users/Sijie/GoogleDrive/PERMIT/3-Simulation Code/packet size obfuscation'+
      '/python_code/run_on_cluster/DiscreteTime/')
sys.path.append(wd)

#%%
import numpy as np
import helper as hp
import NetComponents as nc
from scipy.io import savemat, loadmat

import pdb

#%%
def main():
    ns = hp.NetSetting(var_out_interval=False) 
    
    prefix = 'DPS_VSCI'
    subfolder = 'DPS_vs_VSCI'
    fn = ns.get_fn(subfolder=subfolder, prefix=prefix) # get filename of optim res
    # local
#    subfolder = 'DPS_vs_VSCI'
#    fp = ns.get_local_fp(subfolder=subfolder)
    # cluster
#    fp = '/home/sx37/DT/DPS/data/IoT'
    
    pdb.set_trace()
    mat_res = loadmat(fn)
    
    Q_opt_all = mat_res['Q_opt_all']
    
    sdist, adist = hp.sz_ia_generator(ns=ns)
    
    delay_all = np.ones((ns.n_eps, ns.n_rho))*np.nan
    rho_actual_all = np.ones((ns.n_eps, ns.n_rho))*np.nan
    qsz_actual_all = np.ones((ns.n_eps, ns.n_rho))*np.nan
    dep_rate_all = np.ones((ns.n_eps, ns.n_rho))*np.nan
    
    for i_eps in range(ns.n_eps):
        for i_rho in range(ns.n_rho):
            Q_opt = Q_opt_all[:,:,i_eps,i_rho]
            queue_stats = nc._obj_fun(x=Q_opt, ns=ns, sdist=sdist, adist=adist,
                                      shaper='DPS', return_qsz=True, return_byte_oh=True,
                                      return_dep_rate=True, return_delay_cdf=False)
            
#            pdb.set_trace()
            
            delay = queue_stats['delay']
            rho_actual = 1/(queue_stats['byte_oh']+1)
            qsz_actual = queue_stats['queue_sz']
            dep_rate = queue_stats['dep_rate']
            
            delay_all[i_eps, i_rho] = delay
            rho_actual_all[i_eps, i_rho] = rho_actual
            qsz_actual_all[i_eps, i_rho] = qsz_actual
            dep_rate_all[i_eps, i_rho] = dep_rate
    
    mat_res['delay_all_DPS'] = delay_all
    mat_res['rho_actual_all_DPS'] = rho_actual_all
    mat_res['qsz_actual_all_DPS'] = qsz_actual_all
    mat_res['dep_rate_all_DPS'] = dep_rate_all
    
    pdb.set_trace()
    rn = ns.get_rn(prefix)
    rp = '/home/sx37/DT/DPS/res/'
    savemat(rp+rn, mat_res)

if __name__ == '__main__':
    main()