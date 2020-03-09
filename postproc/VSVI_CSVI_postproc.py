#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 22:22:15 2020

@author: sijie
"""

#%% Import Packages
import sys
wd = '/Users/Sijie/GoogleDrive/PERMIT/3-Simulation Code/packet size obfuscation/python_code/run_on_cluster/DiscreteTime/'
sys.path.append(wd)

import numpy as np
import functools
import helper as hp
import NetComponents as nc

import pdb
from scipy.io import savemat, loadmat

#%% Network and Simulation Setup
def main():
    ns = hp.NetSetting(var_out_interval=True) 
    
    prefix = 'VSVI_CSVI'
    subfolder = 'DPS_vs_VSVI'
    fn = ns.get_fn(subfolder=subfolder, prefix=prefix) # get filename of optim res
    # local
#    subfolder = 'DPS_vs_VSVI'
#    fp = ns.get_local_fp(subfolder=subfolder)
    # cluster
#    fp = '/home/sx37/DT/DPS/data/'
    
    mat_res = loadmat(fn)
    
    mu_opt_all = mat_res['mu_opt_all']
    exp_pkt_sz_all = mat_res['exp_pkt_sz_all'].squeeze()
    
    sdist, adist = hp.sz_ia_generator(ns=ns)
    
    obj_fun = functools.partial(nc._obj_fun, sdist=sdist, adist=adist, 
                                return_qsz=True, return_qsz_depart=True, 
                                return_rho_actual=True, return_byte_oh=True, 
                                return_dep_rate=True)
    
    # no need to consider 0 packet size
    x = np.ones((ns.N,1)) # a special case of Q_opt with an all-1 column
    
    delay_all_VSVI = []
    rho_actual_all_VSVI = []
    qsz_actual_all_VSVI = []
    qsz_depart_actual_all_VSVI = []
    dep_rate_all_VSVI = []
    
    delay_all_CSVI = []
    rho_actual_all_CSVI = []
    qsz_actual_all_CSVI = []
    qsz_depart_actual_all_CSVI = []
    dep_rate_all_CSVI = []
    
    for i_rho in range(ns.n_rho):
#        
#        pdb.set_trace()
        
        # Calculate stats for VSVI
        mu_opt = mu_opt_all[i_rho,:]
        
        # VSVI shaper behaves like DPS: with q being a rank-1 matrix.
        q_vsvi = x*mu_opt
        q_vsvi[0,:] = np.r_[1, [0]*(ns.hN-1)]
        queue_stats_VSVI = obj_fun(x=q_vsvi, ns=ns, shaper='DPS')
        
        delay = queue_stats_VSVI['delay']
        delay_all_VSVI.append(delay)
        
        qsz = queue_stats_VSVI['queue_sz']
        qsz_actual_all_VSVI.append(qsz)
        qsz_depart = queue_stats_VSVI['queue_sz_depart']
        qsz_depart_actual_all_VSVI.append(qsz_depart)
        
        rho_actual = queue_stats_VSVI['rho_actual']
        rho_actual_all_VSVI.append( rho_actual )
        
        dep_rate = queue_stats_VSVI['dep_rate']
        dep_rate_all_VSVI.append(dep_rate)
        
#        pdb.set_trace()
        
        # Calculate stats for CSVI
        exp_pkt_sz = exp_pkt_sz_all[i_rho]
        # change the out_sizes to a constant size
        out_sizes = ns.out_sizes
        hN = ns.hN
        ns.out_sizes = [exp_pkt_sz]
        ns.hN = 1
        # For CSVI, the output is deterministic, no need to set n_out > 1
        ns.n_out = 1
        
        queue_stats_CSVI = obj_fun(x=x, ns=ns, shaper='CSVI')
        # need to change the out_sizes back since 
        # both VSVI & CSVI are in the same loop
        ns.out_sizes = out_sizes 
        ns.hN = hN
        
        delay = queue_stats_CSVI['delay']
        delay_all_CSVI.append(delay)
        
        qsz = queue_stats_CSVI['queue_sz']
        qsz_actual_all_CSVI.append(qsz)
        qsz_depart = queue_stats_CSVI['queue_sz_depart']
        qsz_depart_actual_all_CSVI.append(qsz_depart)
        
        rho_actual = queue_stats_CSVI['rho_actual']
        rho_actual_all_CSVI.append( rho_actual )
        
        dep_rate = queue_stats_CSVI['dep_rate']
        dep_rate_all_CSVI.append(dep_rate)
    
#    pdb.set_trace()
    
    mat_res['delay_all_VSVI'] = delay_all_VSVI
    mat_res['rho_actual_all_VSVI'] = rho_actual_all_VSVI
    mat_res['qsz_actual_all_VSVI'] = qsz_actual_all_VSVI
    mat_res['qsz_depart_actual_all_VSVI'] = qsz_depart_actual_all_VSVI
    mat_res['dep_rate_all_VSVI'] = dep_rate_all_VSVI
    
    mat_res['delay_all_CSVI'] = delay_all_CSVI
    mat_res['rho_actual_all_CSVI'] = rho_actual_all_CSVI
    mat_res['qsz_actual_all_CSVI'] = qsz_actual_all_CSVI
    mat_res['qsz_depart_actual_all_CSVI'] = qsz_depart_actual_all_CSVI
    mat_res['dep_rate_all_CSVI'] = dep_rate_all_CSVI
#    pdb.set_trace()
    
    rn = ns.get_rn(prefix)
    rp = '/home/sx37/DT/DPS/res/'
    savemat(rp+rn, mat_res)

if __name__ == "__main__":
    main()