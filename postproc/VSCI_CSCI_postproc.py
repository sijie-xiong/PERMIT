#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 12:45:39 2020

@author: Sijie
"""

#%% Import Packages
import sys
wd = '/Users/Sijie/GoogleDrive/PERMIT/3-Simulation Code/packet size obfuscation/python_code/run_on_cluster/DiscreteTime/'
sys.path.append(wd)

import numpy as np
import functools
import helper as hp
import NetComponents as nc

#import pdb
from scipy.io import savemat, loadmat

#%% Network and Simulation Setup
def main():
    ns = hp.NetSetting(var_out_interval=False) 
    
    prefix = 'VSCI_CSCI'
    subfolder = 'DPS_vs_VSCI'
    fn = ns.get_fn(subfolder=subfolder, prefix=prefix) # get filename of optim res
    # local
#    subfolder = 'DPS_vs_VSCI'
#    fp = ns.get_local_fp(subfolder=subfolder)
    # cluster
#    fp = '/home/sx37/DT/DPS/data/'
    
    mat_res = loadmat(fn)
    
    mu_opt_all = mat_res['mu_opt_all']
    exp_pkt_sz_all = mat_res['exp_pkt_sz_all'].squeeze()
    
    sdist, adist = hp.sz_ia_generator(ns=ns)
    
    obj_fun = functools.partial(nc._obj_fun, sdist=sdist, adist=adist, shaper='DPS')
    
    # no need to consider 0 packet size
    x = np.ones((ns.N,1)) # a special case of Q_opt with an all-1 column
    
    delay_all_VSCI = []
    rho_actual_all_VSCI = []
    qsz_actual_all_VSCI = []
    dep_rate_all_VSCI = []
    
    delay_all_CSCI = []
    rho_actual_all_CSCI = []
    qsz_actual_all_CSCI = []
    dep_rate_all_CSCI = []
    
    for i_rho in range(ns.n_rho):
#        
#        pdb.set_trace()
        
        # Calculate stats for VSCI
        mu_opt = mu_opt_all[i_rho,:]
        
        # VSCI shaper behaves like DPS: with q being a rank-1 matrix.
        queue_stats_VSCI = obj_fun(x=x*mu_opt, ns=ns)
        
        delay = queue_stats_VSCI['delay']
        delay_all_VSCI.append(delay)
        
        qsz = queue_stats_VSCI['queue_sz']
        qsz_actual_all_VSCI.append(qsz)
        
        rho_actual = queue_stats_VSCI['rho_actual']
        rho_actual_all_VSCI.append( rho_actual )
        
        dep_rate = queue_stats_VSCI['dep_rate']
        dep_rate_all_VSCI.append(dep_rate)
        
#        pdb.set_trace()
        
        # Calculate stats for CSCI
        exp_pkt_sz = exp_pkt_sz_all[i_rho]
        # change the out_sizes to a constant size
        out_sizes = ns.out_sizes
        hN = ns.hN
        ns.out_sizes = [exp_pkt_sz]
        ns.hN = 1
        # For CSCI, the output is deterministic, no need to set n_out > 1
        ns.n_out = 1 
        
        queue_stats_CSCI = obj_fun(x=x, ns=ns)
        # need to change the out_sizes back since 
        # both VSCI & CSCI are in the same loop
        ns.out_sizes = out_sizes 
        ns.hN = hN
        
        delay = queue_stats_CSCI['delay']
        delay_all_CSCI.append(delay)
        
        queue_sz = queue_stats_CSCI['queue_sz']
        qsz_actual_all_CSCI.append(queue_sz)
        
        rho_actual = queue_stats_CSCI['rho_actual']
        rho_actual_all_CSCI.append( rho_actual )
        
        dep_rate = queue_stats_CSCI['dep_rate']
        dep_rate_all_CSCI.append(dep_rate)
    
#    pdb.set_trace()
    
    mat_res['delay_all_VSCI'] = delay_all_VSCI
    mat_res['rho_actual_all_VSCI'] = rho_actual_all_VSCI
    mat_res['qsz_actual_all_VSCI'] = qsz_actual_all_VSCI
    mat_res['dep_rate_all_VSCI'] = dep_rate_all_VSCI
    
    mat_res['delay_all_CSCI'] = delay_all_CSCI
    mat_res['rho_actual_all_CSCI'] = rho_actual_all_CSCI
    mat_res['qsz_actual_all_CSCI'] = qsz_actual_all_CSCI
    mat_res['dep_rate_all_CSCI'] = dep_rate_all_CSCI
#    pdb.set_trace()
    
    rn = ns.get_rn(prefix)
    rp = '/home/sx37/DT/DPS/res/'
    savemat(rp+rn, mat_res)

if __name__ == "__main__":
    main()