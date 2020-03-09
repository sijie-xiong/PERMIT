#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 18:08:33 2019

@author: Sijie
"""

"""
    A bit more detailed set of components to use in packet switching
    queueing experiments.
    Copyright 2014 Greg M. Bernstein
    Released under the MIT license
"""
import numpy as np
import random
import copy
from heapq import heappush, heappop

import simpy
import NetComponents as nc

import pdb

#%% Delay objective function evaluated at each x
def _obj_fun(x=None, ns=None, sdist=None, adist=None, shaper='DPS', 
             return_delay=True, return_delay_cdf=False, 
             return_qsz=True, return_qsz_depart=False, return_qsz_cdf=False,
             return_byte_oh=True, return_rho_actual=True, 
             return_dep_rate=True):
    
#    pdb.set_trace()
    
    in_sizes, out_sizes = ns.in_sizes, ns.out_sizes
    N, hN = len(in_sizes), len(out_sizes)
    
#    pdb.set_trace()
    
    Q_opt = x
    if hN > 1 and x.shape != (N, hN): # if not constant-size departures
        Q_opt = x.reshape(newshape=(N, hN), order='F')
    Q_opt_dict = {}
    for i in range(N):
        Q_opt_dict[in_sizes[i]] = Q_opt[i,:]
    
    def monitor_interval():
        return 1
    
    delay_all, qsz_actual_all, qsz_depart_actual_all, byte_oh_all, rho_actual_all = [],[],[],[],[]
    dep_rate_all, perpkt_delay_all, qsz_all = [], [], []
    for i_rep in range(ns.n_reps):
        
        pkt_list = [sdist() for i in range(int(ns.T/ns.in_tau))]
        
        # Instantiate Network Components
        env = simpy.Environment()
        pg = nc.PacketGenerator(env, adist=adist, pkt_list=pkt_list)
        
        snoop = nc.SnoopSplitter(n_out=ns.n_out)
        pg.out = snoop
        
#        pdb.set_trace()
        
        ts_all = []
        ts_monitor_all = []
        for i_out in range(ns.n_out):
            if shaper == 'DPS':
                ts = nc.DPS(env, Q_opt_dict=Q_opt_dict, out_sizes=out_sizes)
            elif shaper == 'CSCI':
                ts = nc.CSCI(env, out_size=out_sizes[0], out_tau=ns.out_tau)
            elif shaper == 'CSVI':
                ts = nc.CSVI(env, out_size=out_sizes[0])
            else:
                raise ValueError("'shaper' must be one of 'CSCI', 'CSVI' and 'DPS'!")
            ts_monitor = nc.PortMonitor(env, ts, monitor_interval)
            ts_all.append(ts)
            ts_monitor_all.append(ts_monitor)
            setattr(snoop, 'out'+str(i_out), ts)

        env.run(ns.T+ns.in_tau/2)
        
#        pdb.set_trace()
        
        delays, queue_szs, queue_szs_depart, byte_ohs, rho_actuals, dep_rates = [], [], [], [], [], []
        for i_out in range(ns.n_out):
            
            ts = ts_all[i_out]
            ts_monitor = ts_monitor_all[i_out]
            delays.append( ts.waits / ts.packets_rec ) # delay per pkt
            queue_szs.append( np.mean(ts_monitor.sizes) )
            queue_szs_depart.append( np.mean(ts.queue_sizes) )
            byte_oh = ts.byte_overhead / ts.byte_rec 
            byte_ohs.append( byte_oh ) 
            rho_actuals.append( ts.byte_rec / ts.byte_dep )
            dep_rates.append( ts.packets_dep / ns.T )
            
        if return_delay_cdf:
            perpkt_delay_all = np.r_[perpkt_delay_all, ts.perpkt_delays]
        if return_qsz_cdf:
            qsz_all = np.r_[qsz_all, ts_monitor.sizes]
        
        delay_all.append( np.mean(delays) )
        qsz_actual_all.append( np.mean(queue_szs) )
        qsz_depart_actual_all.append( np.mean(queue_szs_depart) )
        byte_oh_all.append( np.mean(byte_ohs) )
        rho_actual_all.append( np.mean(rho_actuals) )
        dep_rate_all.append( np.mean(dep_rates) )
    
#    pdb.set_trace()
    
    # avg queue stats for all epsilons
    avg_delay = np.mean(delay_all)
    avg_queue_sz = np.mean(qsz_actual_all)
    avg_queue_sz_depart = np.mean(qsz_depart_actual_all)
    avg_byte_oh = np.mean(byte_oh_all)
    avg_rho_actual = np.mean(rho_actual_all)
    avg_dep_rate = np.mean(dep_rate_all)
    
    queue_stats = {}
    if return_delay:
        queue_stats['delay'] = avg_delay
    if return_delay_cdf:
        sorted_perpkt_delays = np.sort(perpkt_delay_all)
        cdf = np.arange(1, len(sorted_perpkt_delays)+1)/float(len(sorted_perpkt_delays))
        queue_stats['delay_cdf'] = (sorted_perpkt_delays, cdf)
    if return_qsz:
        queue_stats['queue_sz'] = avg_queue_sz
    if return_qsz_depart:
        queue_stats['queue_sz_depart'] = avg_queue_sz_depart
    if return_qsz_cdf:
        sorted_queue_sizes = np.sort(qsz_all)
        cdf = np.arange(1, len(sorted_queue_sizes)+1)/float(len(sorted_queue_sizes))
        queue_stats['qsz_cdf'] = (sorted_queue_sizes, cdf)
    if return_byte_oh:
        queue_stats['byte_oh'] = avg_byte_oh
    if return_rho_actual:
        queue_stats['rho_actual'] = avg_rho_actual
    if return_dep_rate:
        queue_stats['dep_rate'] = avg_dep_rate
    return queue_stats

#%%
class Shaper(object):
    
    def __init__(self, env):
        self.env = env
        self.store = []
#        self.store = simpy.PriorityStore(self.env)
        
        self.packets_rec = 0        # pkts recorded
        self.byte_rec = 0           # bytes recorded
        self.packets_dep = 0        # pkts departed
        self.byte_dep = 0           # bytes departed
        self.byte_size = 0          # queue size in bytes
        self.waits = 0              # tot waiting time of all pkts
        self.perpkt_delays = []     # delays of each pkt
        self.queue_sizes = []
        self.byte_overhead = 0      # dummy bytes
        
        self.action = self.env.process(self.run())
        self.out = None
        self.out_time = 0
    
    def get_dep_size(self):
#        pdb.set_trace()
        if (hasattr(self, 'out_sizes') and 
            hasattr(self, 'Q_opt_dict') and 
            hasattr(self, 'last_sz')):
            sz = random.choices(self.out_sizes, weights=self.Q_opt_dict[self.last_sz])[0]
        elif hasattr(self, 'out_size'):
            sz = self.out_size
        else:
            raise AttributeError('DPS mechanism needs out_sizes, Q_opt_dict and last_sz.\n'+
                                 'CSCI mechanism needs out_size.\n'+
                                 'CSVI mechanism needs out_size.')
        return sz
        
    def run(self):
#        yield self.env.timeout(0)
        while True:
            # DPS and CSVI shapers have arrival attribute as a SimPy event().
            # - DPS arrival gets reset in every time slot of length in_tau;
            # - CSVI arrival gets reset whenever there's a non-zero arrival.
            if hasattr(self, 'arrival'):
                yield self.arrival
            # CSCI shaper doesn't need/have arrival attribute since pkts depart 
            # at regular intervals of length out_tau.
            elif hasattr(self, 'out_tau'):
                if self.out_tau is not None:
                    yield self.env.timeout(self.out_tau)
                else:
                    raise ValueError("'out_tau' attribute is None, not specified for CSCI!")
            else:
                raise AttributeError("To issue a departure, this shaper must "+
                                     "have attribute of either 'arrival' (DPS/CSVI) "+ 
                                     "or 'out_tau' (CSCI)!")
            
            now = self.env.now
            
#            time_diff = now - self.out_time
#            if time_diff > 1.5:
#                print(self.byte_dep, time_diff, self.byte_dep/now)
#                pdb.set_trace()
#            self.out_time = now
            
            sz = self.get_dep_size()
            if sz > 0:
                self.packets_dep += 1
            self.byte_dep += sz
            
            msg_list = []
            msg = None
            tmp_size = 0
            
#            avg_delay_1 = 0 if self.packets_rec==0 else self.waits/self.packets_rec
#            avg_delay_2 = np.mean(np.r_[self.perpkt_delays, 
#                                        [0]*(self.packets_rec-len(self.perpkt_delays))]) if self.perpkt_delays else 0
#            print('Byte overhead: ', self.byte_overhead)
#            print('Perpkt delays: ', self.perpkt_delays)
#            print('Avg delay 1: ', avg_delay_1)
#            print('Avg delay 2: ', avg_delay_2)
#            
#            print('Time: ', now, ', Out: ', sz)
#            pdb.set_trace()
            
            while tmp_size < sz and len(self.store):
                
#                pdb.set_trace()
                msg = heappop(self.store)
#                msg = yield self.store.get()
                
                if msg.size == 0:
                    continue
                    
                msg_list.append(msg)
        
                tmp_size += msg.size
                self.byte_size -= msg.size
                
                # If the pkt was split,
                if msg.split:
                    # and its remaining size > sz, can split further
                    if msg.size > sz:
                        continue
                    # and its remaining size <= sz, count its delay
                    else:
                        delay = now - msg.time
                        self.waits += delay
                        self.perpkt_delays.append(delay)
                # If the pkt wasn't split, temporarily count its delay,
                # this accounts for the case when sz > a whole pkt size.
                else:
                    delay = now - msg.time
                    self.waits += delay
                    self.perpkt_delays.append(delay)
#                print(now - msg.time, self.waits)
#                pdb.set_trace()
            
            # If the last dequeue packet was split
            rem_size = tmp_size - sz
            if rem_size > 0:
                rem_pkt = Packet(time=msg.time, size=rem_size, id=msg.id, 
                                 split=True, typ=msg.type, 
                                 flow_id=msg.flow_id, src=msg.src, dst=msg.dst)
#                self.store.put(rem_pkt)
                heappush(self.store, rem_pkt)
                # Reset the queue length
                self.byte_size += rem_size
                # Reset the total delay when only part of this packet left
                if msg.split and msg.size > sz:
                    continue
                else:
                    self.waits -= now - msg.time
                    self.perpkt_delays.pop()
#                print(now - msg.time, self.waits)
#                pdb.set_trace()
            # Else dummy packet sent
            else:
                self.byte_overhead += (-rem_size)
            
            # This only append queue size immediately after departure times.
            self.queue_sizes.append(self.byte_size) 
    
#%%
class DPS(Shaper):
    
    def __init__(self, env, Q_opt_dict=None, out_sizes=None):
        Shaper.__init__(self, env)
        
        # Member values unique to this shaper
        self.last_sz = 0            # size of the last arriving packet
        self.Q_opt_dict = Q_opt_dict
        self.out_sizes = out_sizes
        self.arrival = self.env.event()
            
    def put(self, pkt):
        sz = pkt.size
        self.last_sz = sz
        
#        print('Time: ', self.env.now, ', In: ', sz)
#        pdb.set_trace()
        
        if sz > 0:
            self.packets_rec += 1
            self.byte_rec += sz
            self.byte_size += sz
            
            heappush(self.store, pkt)
        
        self.arrival.succeed()
        self.arrival = self.env.event()
        
        return 
    
 #%%
class CSVI(Shaper):
    
    def __init__(self, env, out_size=None):
        Shaper.__init__(self, env)
        
        # Member values unique to this shaper
        self.out_size = out_size
        self.arrival = self.env.event()

    def put(self, pkt):
        sz = pkt.size
        
#        print('Time: ', self.env.now, ', In: ', sz)
#        pdb.set_trace()
        
        if sz > 0:
            self.packets_rec += 1
            self.byte_rec += sz
            self.byte_size += sz
            
            heappush(self.store, pkt)
            
            self.arrival.succeed()
            self.arrival = self.env.event()
            
        return 
    
#%%
class CSCI(Shaper):
    
    def __init__(self, env, out_size=None, out_tau=None):
        Shaper.__init__(self, env)
        
        # Member values unique to this shaper
        self.out_size = out_size
        self.out_tau = out_tau
    
    def put(self, pkt):
        sz = pkt.size
        
#        print('Time: ', self.env.now, ', In: ', sz)
#        pdb.set_trace()
        
        if sz > 0:
            self.packets_rec += 1
            self.byte_rec += sz
            self.byte_size += sz
            
            heappush(self.store, pkt)
        return 
              
#%%
class Packet(object):
    """ A very simple class that represents a packet.
        This packet will run through a queue at a switch output port.
        We use a float to represent the size of the packet in bytes so that
        we can compare to ideal M/M/1 queues.

        Parameters
        ----------
        time : float
            the time the packet arrives at the output queue.
        size : float
            the size of the packet in bytes
        id : int
            an identifier for the packet
        src, dst : int
            identifiers for source and destination
        flow_id : int
            small integer that can be used to identify a flow
    """
    def __init__(self, time, size, id, split=False, typ=None,
                 flow_id=0, src="a", dst="z"):
        self.time = time
        self.size = size
        self.id = id
        self.type = typ
        self.split = split
        self.flow_id = flow_id
        self.src = src
        self.dst = dst
       
    def __lt__(self, other):
        if self.time == other.time:
            return self.type == 'data'
        else:
            return self.time < other.time
            
    def __repr__(self):
        return "id: {}, src: {}, time: {}, size: {}, split: {}".\
            format(self.id, self.src, self.time, self.size, self.split)

#%%
class PacketGenerator(object):
    """ Generates packets with given inter-arrival time distribution.
        Set the "out" member variable to the entity to receive the packet.

        Parameters
        ----------
        env : simpy.Environment
            the simulation environment
        adist : function
            a no parameter function that returns the successive inter-arrival times of the packets
        sdist : function
            a no parameter function that returns the successive sizes of the packets
        initial_delay : number
            Starts generation after an initial delay. Default = 0.
        duration : number
            Stops generation after the duration. Default = inf.


    """
    def __init__(self, env, id=0, adist=None, sdist=None, pkt_list=None,
                 typ=None, initial_delay=0, duration=float("inf"), flow_id=0):
        self.pkt_list = iter(pkt_list)
        self.id = id
        self.env = env
        self.adist = adist
        self.sdist = sdist
        self.type = typ
        self.initial_delay = initial_delay
        self.duration = duration
        self.out = None
        self.packets_sent = 0
        self.action = env.process(self.run())
        self.flow_id = flow_id

    def run(self):
        """The generator function used in simulations.
        """
#        yield self.env.timeout(self.initial_delay)
        while self.env.now < self.initial_delay + self.duration:
            # wait for next transmission
            yield self.env.timeout(self.adist())
            if self.pkt_list:
                sz = next(self.pkt_list)
            else:
                sz = self.sdist()
            if sz > 0:
                self.packets_sent += 1
            p = Packet(self.env.now, sz, self.packets_sent,
                       typ=self.type, src=self.id, flow_id=self.flow_id)
            self.out.put(p)

#%%   
class SnoopSplitter(object):
    """ A snoop port like splitter. Sends the original packet out port 1
        and sends a copy of the packet out port 2.

        You need to set the values of out1 and out2.
    """
    def __init__(self, n_out=None):
        self.n_out = n_out
        for i in range(n_out):
            setattr(self, 'out'+str(i), None)

    def put(self, pkt):
        for i in range(self.n_out):
            if getattr(self, 'out'+str(i)):
                pkt_copy = copy.copy(pkt)
                getattr(self, 'out'+str(i)).put(pkt_copy)
            
#%%
class PacketSink(object):
    """ Receives packets and collects delay information into the
        waits list. You can then use this list to look at delay statistics.

        Parameters
        ----------
        env : simpy.Environment
            the simulation environment
        debug : boolean
            if true then the contents of each packet will be printed as it is received.
        rec_arrivals : boolean
            if true then arrivals will be recorded
        absolute_arrivals : boolean
            if true absolute arrival times will be recorded, otherwise the time between consecutive arrivals
            is recorded.
        rec_waits : boolean
            if true waiting time experienced by each packet is recorded
        selector: a function that takes a packet and returns a boolean
            used for selective statistics. Default none.

    """
    def __init__(self, env, rec_sizes=False, rec_arrivals=False, absolute_arrivals=False, 
                 rec_waits=True, debug=False, selector=None):
        self.store = simpy.Store(env)
        self.env = env
        self.rec_sizes = rec_sizes
        self.rec_waits = rec_waits
        self.rec_arrivals = rec_arrivals
        self.absolute_arrivals = absolute_arrivals
        self.sizes = []
        self.waits = []
        self.arrivals = []
        self.debug = debug
        self.packets_rec = 0
        self.bytes_rec = 0
        self.selector = selector
        self.last_arrival = 0.0

    def put(self, pkt):
        if not self.selector or self.selector(pkt):
            now = self.env.now
            if self.rec_sizes:
                self.sizes.append(pkt.size)
            if self.rec_waits:
                self.waits.append(self.env.now - pkt.time)
            if self.rec_arrivals:
                if self.absolute_arrivals:
                    self.arrivals.append(now)
                else:
                    self.arrivals.append(now - self.last_arrival)
                self.last_arrival = now
            self.packets_rec += 1
            self.bytes_rec += pkt.size
            if self.debug:
                print(pkt)

#%%
class PortMonitor(object):
    """ A monitor for an SwitchPort. Looks at the number of items in the SwitchPort
        in service + in the queue and records that info in the sizes[] list. The
        monitor looks at the port at time intervals given by the distribution dist.

        Parameters
        ----------
        env : simpy.Environment
            the simulation environment
        port : SwitchPort
            the switch port object to be monitored.
        dist : function
            a no parameter function that returns the successive inter-arrival times of the
            packets
    """
    def __init__(self, env, port, dist, count_bytes=True):
        self.port = port
        self.env = env
        self.dist = dist
        self.count_bytes = count_bytes
        self.sizes = []
        self.lengths = []
        self.action = env.process(self.run())

    def run(self):
        yield self.env.timeout(0.5)
        while True:
            yield self.env.timeout(self.dist())
            sz = self.port.byte_size
            num = len(self.port.store)
            # If record queue size in bytes
            if self.count_bytes:
                self.sizes.append(sz)
            # Otherwise record queue length in num of pkts
            else:
                self.lengths.append(num)