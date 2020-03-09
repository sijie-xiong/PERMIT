#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
Configure network and shaper parameters.

Functions
---------
- sz_ia_generator : callable functions for packet sizes and interarrival times.
- get_net_setting : returns (in_sizes, in_rates, out_sizes) of shaper.
- eps_actual : actual privacy level of an optimization-derived shaper Q_opt.
- FindOptQ : finds byte-overhead-minimal shaper Q.

Created on Tue Apr 30 13:04:31 201
@author: Sijie
"""

# Import Packages
import numpy as np
import random
from math import factorial
import functools
from scipy.optimize import linprog
from scipy.io import loadmat

import pdb
        
#%%
class SimRes(object):
    def __init__(self, status=None, Q_opt=None, eps=None, eps_actual=None, 
                 rho=None, rho_est=None):
        self.status = status
        self.Q_opt = Q_opt
        self.eps = eps
        self.eps_actual = eps_actual
        self.rho = rho
        self.rho_est = rho_est
    
    def set_qstats(self, rho_actual=None, 
                   qsz=None, qsz_actual=None, 
                   dep_rate=None, dep_rate_actual=None, 
                   delay=None, delay_nmld=None, delay_cdf=None):
        self.rho_actual = rho_actual
        self.qsz = qsz
        self.qsz_actual = qsz_actual
        self.dep_rate = dep_rate
        self.dep_rate_actual = dep_rate_actual
        self.delay = delay
        self.delay_nmld = delay_nmld
        self.delay_cdf = delay_cdf
        
    def show(self):
        dec = 3
        print('-----------------------------------')
        print('Optimization Status: ', self.status)
        print('Q_opt = \n', np.round(self.Q_opt, dec))
        print('eps = ', self.eps, 
              ' eps_actual = ', np.round(self.eps_actual, dec))
        print('rho = ', np.round(self.rho, dec), 
              ' rho_est = ', np.round(self.rho_est, dec))
        if hasattr(self, 'rho_actual'):
            print('rho_actual = ', np.round(self.rho_actual, dec))
        if hasattr(self, 'qsz') and hasattr(self, 'qsz_actual'):
            print('qsz = ', np.round(self.qsz, dec), 
                  ' qsz_actual = ', np.round(self.qsz_actual, dec))
        if hasattr(self, 'dep_rate') and hasattr(self, 'dep_rate_actual'):
            print('dep_rate = ', np.round(self.dep_rate, dec),
                  ' dep_rate_actual = ', np.round(self.dep_rate_actual, dec))
        if hasattr(self, 'delay') and hasattr(self, 'delay_nmld'):
            print('delay = ', np.round(self.delay, dec), 
                  ' delay_nmld = ', np.round(self.delay_nmld, dec))
        print('-----------------------------------')
        
#%%
def _Geometric(sizes, p=0.5):
    n_sizes = len(sizes)
    k = np.arange(1, n_sizes+1)
    pdf = (1-p)**(k-1) * p
    pdf /= np.sum(pdf)
    return pdf

def _Poisson(sizes, r=1):
    n_sizes = len(sizes)
    k = np.arange(1, n_sizes+1)
    pdf = np.divide(r**k * np.exp(-r), [factorial(i) for i in k])
    pdf /= np.sum(pdf)
    return pdf

# Zipf distribution
def _Zipf(sizes, s=1):
    n_sizes = len(sizes)
    denom = 0
    for j in range(n_sizes):
        denom += (1/(j+1))**s
    pdf = 1 / np.arange(1, n_sizes+1)**s / denom
    return pdf

class NetSetting():

    def __init__(self, 
                 synthetic=False, device='SleepCamMerged', 
                 where='local', input_typ='iid',
                 n_dev=2, n_eps=5, n_rho=6, 
                 pmf='Zipf', scale=0.01, 
                 T=100000, in_tau=1, out_tau=None,
                 n_reps=10, n_out=5, 
                 same_out_sizes=True, 
                 var_out_interval=False, 
                 scaled_alphabet=False, 
                 increased_alphabet=False,
                 how_to_increase='same_cv'):
        # Set default values
        args = locals()
        for attr, val in list(args.items()):
            if attr in ['self','args']:
                continue
            setattr(self, attr, val)
        
        eps_pool = np.array([0.01, 0.1, 1, 2, 5, 10])
        if n_eps > len(eps_pool):
            raise Exception('Not enough eps values in the pool!')
        self.eps_all = eps_pool[:n_eps]
        
        if not self.synthetic:
            self.set_data()
        self.in_sizes = self.get_in_sizes()
        if scaled_alphabet:
            self.in_sizes *= 3
        self.in_rates = self.get_in_rates()
        self.in_byte_rate = self.get_mu()
        self.coeff_var = self.get_cv()
        self.arrival_rate = 1-self.in_rates[0]
        
        if same_out_sizes:
            self.out_sizes = self.in_sizes
            # min utilization when everything padded to the largest size
            min_rho = self.in_byte_rate / np.max(self.out_sizes)
            if var_out_interval:
                min_rho /= self.arrival_rate
            self.rho_all = np.linspace(0.98, min_rho+0.01, self.n_rho)
        else:
            out_sizes, rho_all = [0], []
            # min utilization when everything padded to the largest size
            min_rho = self.in_byte_rate / np.max(self.in_sizes) / self.arrival_rate
            rho_all_tmp = np.linspace(0.98, min_rho+0.01, n_rho)
            for rho in rho_all_tmp:
                out_size = np.ceil( self.in_byte_rate / rho / self.arrival_rate )
                out_sizes.append(out_size)
                rho_all.append( self.in_byte_rate / out_size / self.arrival_rate )
        
        self.N = len(self.in_sizes) 
        self.hN = len(self.out_sizes)
        
        if increased_alphabet:
            n_dev = 4
            self.get_new_ns(n_dev=n_dev)
    
    def set_data(self):
        data_fn = self.device + '.mat'
        if self.where == 'local':
            path = ('/Users/Sijie/GoogleDrive/PERMIT/3-Simulation Code' +
                    '/packet size obfuscation/python_code/other/SmartHomeIoT/' +
                    'matfiles/Filtered/')
            self.data = loadmat(path+data_fn)
        elif self.where == 'cluster':
            path = '/home/sx37/DT/DPS/data/IoT/'
            self.data = loadmat(path+data_fn)
        else:
            raise ValueError("Choose 'where' from 'local' or 'cluster'!")
        if 'trace' in self.data:
            self.trace = self.data['trace'].squeeze()
            self.trace_iter = iter(self.trace)
        if self.input_typ == 'bursty':
            self.T = len(self.trace)
            self.n_reps = 1
            
    def get_in_sizes(self, n_dev=None):
        if self.synthetic:
            log2_minsz = 5
            if n_dev is not None:
                nn_dev = n_dev
            else:
                nn_dev = self.n_dev
            in_sizes = np.r_[0,2**np.linspace(log2_minsz,log2_minsz+nn_dev-1, nn_dev)]
        else:
            in_sizes = self.data['sizes'].squeeze()
        return in_sizes
        
    def get_in_rates(self, in_sizes=None, scale=None):
        if self.synthetic:
            if (in_sizes is None) or (scale is None):
                in_sizes, scale = self.in_sizes, self.scale
            if self.pmf == 'Zipf':
                in_rates = _Zipf(in_sizes, scale)
            elif self.pmf == 'Poisson':
                in_rates = _Poisson(in_sizes, scale);
            else:
                raise ValueError("Choose 'pmf' from 'Zipf' or 'Poisson'!")
        else:
            in_rates = self.data['pmf'].squeeze()
        return in_rates
        # in_rates = _Zipf(in_sizes, s=5)   # 5, 1, 0.01
        # s = 5,    cv = 5.54  ^
        # s = 1,    cv = 1.21  | lighter traffic
        # s = 0.01, cv = 0.82  |
        # in_rates = _Poisson(in_sizes, r=1)[-n_dev:] 
        # r = 0.05, cv = 6.32
        # r = 0.1,  cv = 4.67
        # r = 1,    cv = 1.34
        # r = 3,    cv = 0.69
        # in_rates = _Geometric(in_sizes, p=0.9)[-n_dev:] 
        # p = 0.99, cv = 10.36
        # p = 0.5,  cv = 1.58
        # p = 0.1,  cv = 1.02
    
    def get_mu(self, sizes=None, rates=None):
        if (sizes is None) or (rates is None):
            sizes, rates = self.in_sizes, self.in_rates
        if np.sum(rates) != 1:
            rates /= np.sum(rates)
        mu = np.dot(sizes, rates)
        return mu
    
    def get_cv(self, sizes=None, rates=None):
        if (sizes is None) or (rates is None):
            sizes, rates = self.in_sizes, self.in_rates
        mu = self.get_mu(sizes=sizes, rates=rates)
        if np.sum(rates) != 1:
            rates /= np.sum(rates)
        cv = np.sqrt(np.dot((sizes - mu)**2, rates)) / mu
        return cv
    
    def get_new_ns(self, n_dev=None):
        n = 5000
        new_in_sizes = self.get_in_sizes(n_dev=n_dev)
        scale_all = np.linspace(0.01,10,n)
        Bin_all = []
        CV_all = []
        for new_scale in scale_all:
            new_in_rates = self.get_in_rates(new_in_sizes, new_scale)
            Bin_all.append(self.get_mu(new_in_sizes, new_in_rates))
            CV_all.append(self.get_cv(new_in_sizes, new_in_rates))
        if self.how_to_increase == 'same_bin':
            closest_idx = np.argmin(np.abs(self.in_byte_rate-Bin_all))
        elif self.how_to_increase == 'same_cv':
            closest_idx = np.argmin(np.abs(self.coeff_var-CV_all))
        else:
            raise ValueError("Choose 'how_to_increase' from 'same_bin' or 'same_cv'!")
        closest_scale = scale_all[closest_idx]
        ns = NetSetting(n_dev=n_dev, n_eps=self.n_eps, n_rho=self.n_rho, 
                        pmf=self.pmf, scale=closest_scale,
                        same_out_sizes=self.same_out_sizes,
                        var_out_interval=self.var_out_interval,
                        scaled_alphabet=self.scaled_alphabet,
                        increased_alphabet=False, # temporarily set as false to avoid infinite loop
                        how_to_increase=self.how_to_increase)
        for attr in ['n_dev','scale','in_sizes','in_rates',
                     'in_byte_rate','coeff_var','arrival_rate',
                     'out_sizes','N','hN']:
            setattr(self, attr, getattr(ns,attr))
        return
    
    # Local path for optimization result
    def get_fp(self, subfolder=None):
        if self.synthetic:
            if self.where == 'local':
                fp = ('/Users/Sijie/GoogleDrive/PERMIT/3-Simulation Code' +
                      '/packet size obfuscation/python_code/run_on_cluster/' +
                      'DiscreteTime/matlab_code/res/linspace_rho_all/%s/') % subfolder
            if self.where == 'cluster':
                fp = '/home/sx37/DT/DPS/data/'
        else:
            if self.where == 'local':
                fp = ('/Users/Sijie/GoogleDrive/PERMIT/3-Simulation Code' +
                      '/packet size obfuscation/python_code/run_on_cluster/' +
                      'DiscreteTime/matlab_code/res/IoT/')
            if self.where == 'cluster':
                fp = '/home/sx37/DT/DPS/data/IoT/'
        return fp
    # File name suffix
    def get_format(self):
        if self.synthetic:
            fmt = (('_%d_dev_%.2f_scale_%2.2f_ibr_%.2f_ar_%.2f_cv_%d_%d.mat') % 
                     tuple(np.r_[self.n_dev, self.scale, self.in_byte_rate, self.arrival_rate, 
                                 self.coeff_var,self.n_eps,self.n_rho]))
        else:
            fmt = (('_%d_szs_%2.2f_ibr_%.2f_ar_%.2f_cv_%d_%d.mat') % 
                     tuple(np.r_[self.N, self.in_byte_rate, self.arrival_rate, 
                                 self.coeff_var,self.n_eps,self.n_rho]))
        return fmt
    def get_fn(self, subfolder=None, prefix=None):
        fp = self.get_fp(subfolder=subfolder)
        fmt = self.get_format()
        fn = prefix+'_opt'+fmt
        if not self.synthetic:
            fn = self.device+'_'+fn
        return fp+fn
    def get_rn(self, prefix=None):
        fmt = self.get_format()
        fn = prefix+'_'+self.input_typ+fmt
        if not self.synthetic:
            fn = self.device+'_'+fn
        return fn
    
    def __repr__(self):
        res = ''
        for attr, val in self.__dict__.items():
            res += attr+': '+str(val)+'\n'
        return res
#        pdb.set_trace()
                                        
#        npts = 9
#        p1_vec = np.linspace(0.1,0.9,npts)
#        p2_vec = np.linspace(0.1,0.9,npts)
#        coeff_var_all = np.ones((npts,npts))*np.nan
#        for ip1 in range(len(p1_vec)):
#            for ip2 in range(len(p2_vec)):
##                pdb.set_trace()
#                p1 = p1_vec[ip1]
#                p2 = p2_vec[ip2]
#                p0 = 1-(p1+p2)
#                if p0 > 0 and p0 < 1:
#                    in_rates = np.r_[p0,p1,p2]
#                    in_byte_rate = np.dot(self.in_sizes, in_rates)
#                    coeff_var = np.sqrt(np.dot((self.in_sizes - in_byte_rate)**2, in_rates)) / in_byte_rate
#                    coeff_var_all[ip1,ip2] = coeff_var
#        print(np.any(coeff_var_all<1))
#        pdb.set_trace()
    
#%%  
# generate variable packet sizes
def _var_sz(sizes=None, rates=None):
    sz = random.choices(sizes, weights=rates)[0]  
    return sz
# generate constant packet size
def _const_sz(size=None):
    return size
def _sz_from_trace(ns=None):
    try:
        return next(ns.trace_iter)
    except:
        ns.trace_iter = iter(ns.trace)
        return _sz_from_trace(ns)

def _const_ia(in_tau=None):
    return in_tau
# generate Exponential inter-arrival times
def _exp_ia(rate=None):
    return np.random.exponential(1/rate)

# Markov chain for generating input traffic
class MC(object):
    def __init__(self, ns=None, trans_mat=None):
        self.current_sz = 0
        self.ns = ns
        self.P = trans_mat
        
        self.P_dict = {}
        for i in range(ns.N):
            self.P_dict[ns.in_sizes[i]] = self.P[i,:]
    
    def transition(self):
        sz = _var_sz(sizes=self.ns.in_sizes, rates=self.P_dict[self.current_sz])
        self.current_sz = sz
        return sz
        
# packet sizes and interarrival times generator
def sz_ia_generator(ns=None):
    """ Callable 'packet-size' and 'interarrival-time' generating functions.
    
    Parameters
    ----------
    sizes : 1-D array-like
        Packet sizes, integers.
    rates : 1-D array-like
        Packet arrival rates, in accordance with sizes.
    in_tau: float
        Constant interval.
    typ: str
        'CSVI' or 'VSCI', type of packet size/interval generator.
    """
    # single pkt generator: 'aggregate' traffic with variable-size/discrete-time pkts
    if ns.input_typ == 'iid':
        sz_gen = functools.partial(_var_sz, sizes=ns.in_sizes, rates=ns.in_rates)
        ia_gen = functools.partial(_const_ia, in_tau=ns.in_tau)
        return sz_gen, ia_gen
    
    if ns.input_typ == 'bursty':
        sz_gen = functools.partial(_sz_from_trace, ns=ns)
        ia_gen = functools.partial(_const_ia, in_tau=ns.in_tau)
        return sz_gen, ia_gen
    
    if ns.input_typ == 'markov':
        # The transition matrix of iid process as a special case of MC, 
        # all rows equal to the iid marginal distribution.
        P_iid = np.kron(np.ones((ns.N,1)), ns.in_rates) 
        # Eigen decomposition of the transition matrix.
        w, v = np.linalg.eig(P_iid.T)
        
        msg1 = "The largest eigen value of transition matrix isn't 1!"
        msg2 = "The stationary distribution doesn't match iid pkt sz distribution!"
        # Assert that the largest eigen value is 1.
        w_max_idx = np.argmax(w)
        assert (np.abs(w[w_max_idx]-1)<1e-10), msg1
        pi = v[:,w_max_idx]/sum(v[:,w_max_idx])
        # Assert that the eigenvector corresponding to the eigenvalue 1 is
        # indeed equal to the iid marginal distribution
        assert (np.sum(np.abs(pi - ns.in_rates))<1e-10), msg2
        
        # Tweak the eigen values a bit to form a new transition matrix, 
        # with which the resulting process isn't iid anymore.
        u = np.r_[1, [0.1]*(ns.N-1)]
        # Reconstruct a transition matrix with updated eigen values.
        P_markov = v.dot(np.diag(u)).dot(np.linalg.inv(v)).T
        w, v = np.linalg.eig(P_markov.T)
        w_max_idx = np.argmax(w)
        assert (np.abs(w[w_max_idx]-1)<1e-10), msg1
        pi = v[:,w_max_idx]/sum(v[:,w_max_idx])
        assert (np.sum(np.abs(pi - ns.in_rates))<1e-10), msg2
        
        mc = MC(ns=ns, trans_mat=P_markov)
        
        sz_gen = mc.transition
        ia_gen = functools.partial(_const_ia, in_tau=ns.in_tau)
        return sz_gen, ia_gen

    # multiple pkt generators: 'individual' traffic, each of constant-size/Exp-interval pkts
    if ns.input_typ == 'ind':
        sdist_all, adist_all = [], []
        for i_dev in range(ns.n_dev):
            sz_gen = functools.partial(_const_sz, size=ns.in_sizes[i_dev])
            ia_gen = functools.partial(_exp_ia, rate=ns.in_rates[i_dev])
            sdist_all.append(sz_gen)
            adist_all.append(ia_gen)
        return sdist_all, adist_all
    
#%%
def main():
    n_eps, n_rho = 5, 6
    # For CSCI, the output is deterministic, no need to set n_out > 1
    ns = NetSetting(n_dev=2, scale=0.01, n_out=1, n_eps=n_eps, n_rho=n_rho, same_out_sizes=True) 
    
    sdist, adist = sz_ia_generator(ns=ns, typ='markov')
    
    pdb.set_trace()
    
    pkt_list = [sdist() for t in range(2*ns.T)]
#    pkt_list = pkt_list[ns.T:]
    
    marginal_freq = np.array([pkt_list.count(in_size) for in_size in ns.in_sizes]) / (2*ns.T)
    print('Freq: ', marginal_freq)
    print('Prior: ', ns.in_rates)
    
    pdb.set_trace()
    
    
if __name__ == "__main__":
    main()
    
#%%
# avg queue size estimate by Wiener Hopf Factorization
def qsz_est_whf(x=None, ns=None, delta=1e-5, alg=2, return_cnt=False):
    Q_opt = x
    
    in_sizes, out_sizes = ns.in_sizes, ns.out_sizes
    N, hN = len(in_sizes), len(out_sizes)
    
    # change into matrix form
    if Q_opt.shape != (N, hN):
        Q_opt = np.reshape( Q_opt, newshape=(N, hN), order='F' )
    
    shaper = OptimalShaper(Q_opt, ns)
    rho_est = shaper.rho_est()
    
#    print(rho_est)
#    pdb.set_trace()
    
    if (np.any(x<0) or np.any(x>1) or 
        np.any(np.sum(Q_opt, axis=1) > 1+1e-6) or
        np.any(np.sum(Q_opt, axis=1) < 1-1e-6) or
        rho_est >= 1 or rho_est < 0):
        return 1e4
    
    ladder_heights = (np.kron(in_sizes.reshape((N,1)), np.ones((1,hN))) - 
                      np.kron(out_sizes.reshape((1,hN)), np.ones((N,1))) )
    
    in_rates = ns.in_rates
    tt = np.kron(in_rates.reshape((N,1)), np.ones((1,hN))) 
    ladder_heights_prob = tt * Q_opt
    
    g = int(np.abs(np.min(ladder_heights)))
    h = int(np.max(ladder_heights))
    
    u_dict = {j:0 for j in range(-g, h+1)}
    
    for i in range(N):
        for j in range(hN):
            ladder = ladder_heights[i,j]
            u_dict[ladder] += ladder_heights_prob[i,j]
    
    a_dict = {i:0 for i in range(1,h+1)}
    cnt = 0
    
    if alg == 1:
        b_dict = {j:0 for j in range(g+1)}
        
        while True:
            cnt += 1
            
            a_dict_new, b_dict_new = {}, {}
            for j in range(g+1):
                b_dict_new[j] = u_dict[-j]
                for i in range(1,h+1):
                    if j+i in b_dict:
                        b_dict_new[j] += a_dict[i] * b_dict[j+i] / (1 - b_dict[0])
                    else:
                        break
#                b_dict_new[j] = min(1, b_dict_new[j])

            max_diff = -np.inf
            for i in range(1,h+1):
                a_dict_new[i] = u_dict[i]
                for j in range(1,g+1):
                    if i+j in a_dict:
                        a_dict_new[i] += a_dict[i+j] * b_dict[j] / (1 - b_dict[0])
                    else:
                        break
                diff = np.abs(a_dict_new[i] - a_dict[i])
                if max_diff < diff:
                    max_diff = diff
                a_dict_new[i] = min(1, a_dict_new[i])
            
            if cnt > 2000:
                print('Utilization: ', rho_est)
                print('Q_opt: \n', Q_opt)
                print('Row sum: ', Q_opt.sum(axis=1))
                raise ValueError('This optimization iterate prevents WHF from converging!')
                
            a_dict, b_dict = a_dict_new, b_dict_new
            
            if max_diff < delta:
                break
        
        S = np.sum(list(b_dict.values())) - b_dict[0]
        for i in a_dict:
            a_dict[i] /= S
    else:
        b_dict = {}
        
        while True:
            cnt += 1
            if cnt > 2000:
                print('Utilization: ', rho_est)
                print('Q_opt: \n', Q_opt)
                print('Row sum: ', Q_opt.sum(axis=1))
                raise ValueError('This optimization iterate prevents WHF from converging!')
            
#            pdb.set_trace()
            
            for j in range(g,-1,-1):
                b_dict[j] = u_dict[-j]
                for i in range(1,h+1):
                    if j+i in b_dict:
                        b_dict[j] += a_dict[i] * b_dict[j+i]
                    else:
                        break
            b_dict[0] = min(1, b_dict[0])
            
#            pdb.set_trace()
            
            if alg == 2:
                S = np.sum(list(b_dict.values())) - b_dict[0] # sum from b_1
            elif alg == 3:
                S = 1-b_dict[0]
                for j in range(1,g+1):
                    S += b_dict[j]
                S = S/2
            
            max_diff = -np.inf
            for i in range(h,0,-1):
                ai = u_dict[i]
                for j in range(1,g+1): # sum from b_1
                    if i+j in a_dict:
                        ai += a_dict[i+j] * b_dict[j]
                    else:
                        break
                ai = ai/S
                
                diff = np.abs(ai-a_dict[i])  
                a_dict[i] = ai
                
                if max_diff < diff:
                    max_diff = diff
            
#            pdb.set_trace()
            
            if max_diff < delta:
                break
    
    w0 = 1-np.sum(list(a_dict.values()))
    queue_sz_est = 0
    for i in range(1,h+1):
        queue_sz_est += i*a_dict[i]
    queue_sz_est /= w0
    
#    print('Total iterations: ', cnt)
#    print('Queue size est.: ', queue_sz_est)
#    pdb.set_trace()
    
    if return_cnt:
        return queue_sz_est, cnt
    return queue_sz_est
    
#%%
def eps_actual(Q_opt=None):
    """ Actual privacy level of an optimization-derived shaper Q_opt.
    """
    Q_opt = np.round(Q_opt, 15)
    N = Q_opt.shape[0]
    
    tt1 = np.dot( np.kron( np.eye(N), np.ones((N,1)) ), Q_opt )
    tt2 = np.dot( np.kron( np.ones((N,1)), np.eye(N) ), Q_opt )

    np.seterr(divide='ignore', invalid='ignore')
    tmp = tt1 / tt2
    np.seterr(divide='warn', invalid='warn')
    
    eps_actual = np.log( np.nanmax(tmp) )
    
#    print('Q_opt: \n', Q_opt)
#    print('eps_actual: ', eps_actual)
#    pdb.set_trace()
    return eps_actual

def dep_rate(Q_opt=None, ns=None):
    marginal = ns.in_rates.dot(Q_opt)
    return 1-marginal[0]
    
class OptimalShaper():
    
    def __init__(self, Q_opt=None, ns=None):
        self.Q_opt = Q_opt
        self.ns = ns
    
    def eps_actual(self):
        return eps_actual(self.Q_opt)
    
    def rho_est(self):
        return self.in_byte_rate()/self.out_byte_rate()
    
    def in_byte_rate(self):
        return np.dot(self.ns.in_rates, self.ns.in_sizes)
    
    def out_byte_rate(self):
        return np.dot(np.dot(self.ns.in_rates, self.Q_opt), self.ns.out_sizes)
    
    def dep_rate(self):
        return dep_rate(self.Q_opt, self.ns)
    
    def qsz_est_whf(self, delta=1e-5, alg=2): 
        return qsz_est_whf(Q_opt=self.Q_opt, ns=self.ns, delta=delta, alg=alg)

#%% Find byte-overhead-optimal Q
def FindOptQ(ns=None, eps=None, rho=None):
    """ Find optimal shaper maatrix Q minimizing byte overhead.
    
    Params
    ------
    in_sizes : 1-D array-like
        Packet sizes.
    in_rates : 1-D array-like
        Corresponding arrival rates of various packet sizes.
    out_sizes : 1-D array-like
        Output packet sizes of the shaper.
    eps : float
        Privacy level.
    rho : float
        Utilization.
    
    Returns
    -------
    opt_shaper : OptimalShaper() object
        
    """
#    pdb.set_trace()
    N, hN = len(ns.in_rates), len(ns.out_sizes)
    
    Byte_ce_mat = np.ravel(
                    np.dot( np.diag( ns.in_rates ),
                            np.kron( np.ones((ns.N, 1)), ns.out_sizes) 
                          ), 
                    order='F' # flatten an array column-wise
                    )
        
    # stability condition
    A_ub1 = -Byte_ce_mat
    b_ub1 = -np.dot(ns.in_rates, ns.in_sizes) / rho
    
    # privacy constraints
    Priv_ci1 = np.kron( 
                np.eye(ns.hN), 
                np.kron( np.eye(ns.N), np.ones((ns.N,1)) ) 
              - np.exp(ns.eps + 1e-3) * np.kron( np.ones((ns.N,1)), np.eye(ns.N) ) 
              )
    A_ub2 = Priv_ci1
    b_ub2 = np.zeros(ns.N**2 * ns.hN)
    # added linear inequality constraint to make privacy level tighter
    Priv_ci2 = np.kron( 
                np.eye(ns.hN), 
                np.kron( np.eye(ns.N), np.ones((ns.N,1)) ) 
              - np.exp(ns.eps - 1e-3) * np.kron( np.ones((ns.N,1)), np.eye(ns.N) ) 
              )
    A_ub3 = -Priv_ci2
    b_ub3 = np.zeros(ns.N**2 * ns.hN)
    
    # stack inequality constraints
    A_ub = np.vstack((A_ub1, A_ub2, A_ub3))
    b_ub = np.r_[b_ub1, b_ub2, b_ub3]
    
    # rows-sum-to-one constraint
    A_eq = np.kron( np.ones((1,hN)), np.eye(N) )
    b_eq = np.ones(N)
    
    bnds = [(0,1)]*(N*hN)
    opts = {'disp': False, 'tol': 1e-6}
    
    res = linprog(Byte_ce_mat, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, 
                  bounds=bnds, method='revised simplex', options=opts)
    
#    pdb.set_trace()
    if res.success:
        Q_opt = np.reshape( res.x, newshape=(N, hN), order='F' )
    else:
        Q_opt = np.nan * np.ones((N, hN))
        
    return OptimalShaper(Q_opt=Q_opt, ns=ns)        
    
def _FindOptQ(ns=None, eps=None, rho=None):
    """ Find optimal shaper Q minimizing byte overhead.
    
    Params
    ------
    in_sizes : 1-D array-like
        Packet sizes.
    in_rates : 1-D array-like
        Corresponding arrival rates of various packet sizes.
    out_sizes : 1-D array-like
        Output packet sizes of the shaper.
    eps : float
        Privacy level.
    rho : float
        Utilization.
    
    Returns
    -------
    opt_shaper : OptimalShaper() object
        
    """
#    pdb.set_trace()
    N, hN = len(ns.in_rates), len(ns.out_sizes)
    
    fun = np.ravel(
            np.dot( np.diag( ns.in_rates ),
                    np.kron( np.ones((N, 1)), ns.out_sizes) 
                  ), 
            order='F' # flatten an array column-wise
            )
    
    # stability condition
    A_ub1 = -fun
    b_ub1 = -np.dot(ns.in_rates, ns.in_sizes) / rho
    
    # privacy constraint
    A_ub2 = np.kron( 
                np.eye(hN), 
                np.kron( np.eye(N), np.ones((N,1)) ) 
              - np.exp(eps) * np.kron( np.ones((N,1)), np.eye(N) ) 
              )
    b_ub2 = np.zeros(N**2 * hN)
    
    # delay constraint (tentative!)
#    n_block = np.min((N, hN))
#    A_ub3 = np.eye(n_block*N)
#    for b in range(n_block): 
#        A_ub3[b*N:(b+1)*N, b*N+b] += -1
#    b_ub3 = np.ones(n_block * N) * -1e-16
#    
#    pdb.set_trace()
#    
#    # stack inequality constraints
#    A_ub = np.vstack((A_ub1, A_ub2, A_ub3))
#    b_ub = np.r_[b_ub1, b_ub2, b_ub3]
    
    # stack inequality constraints
    A_ub = np.vstack((A_ub1, A_ub2))
    b_ub = np.r_[b_ub1, b_ub2]
    
    # rows-sum-to-one constraint
    A_eq = np.kron( np.ones((1,hN)), np.eye(N) )
    b_eq = np.ones(N)
    
    bnds = [(0,1)]*(N*hN)
    opts = {'disp': False, 'tol': 1e-4}
    
    res = linprog(fun, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, 
                  bounds=bnds, method='interior-point', options=opts)
    
    if res.success:
        Q_opt = np.reshape( res.x, newshape=(N, hN), order='F' )
    else:
        Q_opt = np.nan * np.ones((N, hN))
        
    return OptimalShaper(Q_opt=Q_opt, ns=ns)
