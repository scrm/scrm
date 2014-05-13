#!/usr/bin/env python

import numpy as np
import os
import sys
import pylab
from scipy.integrate import simps, trapz

__mydebug__       = False
__fix_ms_seed__   = False
__fix_scrm_seed__ = False

# exact_windows_length [0, 10e2, 10e3, 10e4, -1]

class parameter:
    def __init__( self, nsam = 6, replicate = 100, seqlen = 1e7, rho = 1e-8 , exact_window_length = 0 , divergence = 0, jobs = []): 
        """
        Define simulation parameters
        """
        self.nsam      = nsam
        self.rep       = replicate
        self.seqlen    = seqlen # 10e7
        #self.theta     = 7 * 10e-10 * seqlen * 4 * 10000
        #self.rho       = rho * (self.seqlen-1) * 4 * 10000
        self.rho       = rho
        self.exact_window_length = exact_window_length
        self.divergence = divergence
        self.jobs = jobs
        
        
        delta_points = 20
        big_delta_max = 2e5
        #big_delta_max = 1e4
        small_delta_max = 2e4
        #small_delta_max = big_delta_max
        self.big_delta = np.linspace( 0, int(big_delta_max+1), delta_points )
        #self.big_delta = range( 0, int(big_delta_max+1), 10000 )
        self.small_delta = np.linspace( 0, int(small_delta_max+1), delta_points )

    
    def printing ( self ):
        """
        Check simulation parameters
        """
        print "sample size:       ", self.nsam
        print "replicate:         ",   self.rep
        print "Sequence length:   ", self.seqlen
        print "recombination rate:", self.rho
        print "jobs:", self.jobs
        

def extract_all_info ( job_prefix, num_rep):
    data = []
    for rep in range(1, num_rep+1):
        out = extract_info ( job_prefix, rep )
        data.append( out )
    return data
    
    
def extract_info (job_prefix, ith_rep): 
    prefix = job_prefix + `ith_rep` + "/" + job_prefix + `ith_rep`
    tree_file_name  = prefix + "Trees"
    tree_freq_name  = prefix + "TreeFreq"
    tmrca_name      = prefix + "Tmrca"
    tree_freq = [ float(x) for x in open( tree_freq_name, "r" ) ]    
    tmrca = [ float(x) for x in open( tmrca_name, "r" )]
    timeFile_name = prefix+"time.text"
    runtime = float(open( timeFile_name, "r").readlines()[1].strip("user").strip('\n'))
    return  tree_freq, tmrca, runtime


#def cal_ac_TMRC_star (tree_freq, tmrca, avg_tmrca, delta):

    #seqlen = sum(tree_freq)

    #cumfreq = [0]
    #shifted_cumfreq = [] 
    #for x in tree_freq : 
        #cumfreq.append(cumfreq[-1] + x)
        #shifted_cumfreq.append(cumfreq[-1] - delta)
    #cumfreq.pop(0)

    #n = int(seqlen - delta)
    #ac = 0
    #var = 0
    
    #if __mydebug__:
        #print tree_freq
        #print cumfreq
        #print shifted_cumfreq        
        #print "seqence length = ", seqlen
        #print "i should iterate until ", n, "trees"
    
    #T1_index = 0
    #T2_index = 0 # use while loop to determine the initial T2_index
    #while delta > cumfreq[T2_index]: # need to check again ...
        #T2_index += 1
    #about_to = False
    #for i in range( n ): 
        #if __mydebug__ & ((not i < cumfreq[T1_index]) | (not i < shifted_cumfreq[T2_index])):
            #print "before", i, T1_index,  T2_index
            #print "-------------------------------"
            #about_to = True
        #if i >= cumfreq[T1_index]: T1_index += 1
        #term1 = tmrca[T1_index] - avg_tmrca
        #if i >= shifted_cumfreq[T2_index]: T2_index += 1
        #term2 = tmrca[T2_index] - avg_tmrca
        #if __mydebug__ & about_to:
            #print "after ", i, T1_index,  T2_index
            #about_to = False
        ##print i, T1_index, term1, T2_index, term2
        #ac += term1 * term2
        #var += term1 * term1
    #ac /= float(n)
    #var /= float(n)
    #return ac, var


def cal_ac_TMRC_star_2 (tree_freq, tmrca, avg_tmrca, delta):

    seqlen = int(sum(tree_freq))
    #print seqlen
    cumfreq = [0]
    for x in tree_freq : 
        cumfreq.append(cumfreq[-1] + x)
    cumfreq.pop(0)

    n = int(seqlen - delta)
    ac = 0
    var = 0
    
    T1_index = 0
    T2_index = 0 # use while loop to determine the initial T2_index
    while delta > cumfreq[T2_index]: # need to check again ...
        T2_index += 1
    i = 0
    while i < n:
        distance = min( cumfreq[T1_index] - i, 
                        cumfreq[T2_index] - delta - i )
        distance = min( distance, n - i )
        term1 = tmrca[T1_index] - avg_tmrca
        term2 = tmrca[T2_index] - avg_tmrca
        ac += term1 * term2 * distance
        var += term1 * term1 * distance
        i += distance            
        if i >= cumfreq[T1_index]: T1_index += 1
        if i + delta >= cumfreq[T2_index]: T2_index += 1
    ac /= float(n)
    var /= float(n)
    return ac, var


def process_time (data) :
    run_time_array = []
    for d in data:
        run_time_array.append ( d[2] )
    return np.mean(run_time_array), np.std(run_time_array)
        
def process_data ( data , delta ) :
    # compute the average Tmrca 
    tot_tmrca = 0
    tot_time = 0
    for d in data:
        for i, duration_i in enumerate( d[0] ):   # d[0] : tree_freq, duration of the tree
            tmrca_i = d[1][i]                     # d[1] : tmrca
            tot_time += duration_i                # d[0] : duration of the tree
            tot_tmrca += tmrca_i * duration_i     # JOE changed from tot_tmrca += tmrca_i
    avg_tmrca = tot_tmrca / tot_time
            
    ac_TMRCA  = []    
    ac_mean   = []
    ac_std    = []
    
    for delta_i in delta:
        print "processing delta:", delta_i        
        cum_ac_TMRCA = 0
        cum_var_TMRCA = 0
    
        ratio = []
        for d in data:
            ac, var = cal_ac_TMRC_star_2( d[0], d[1], avg_tmrca, delta_i )
            cum_ac_TMRCA += ac
            cum_var_TMRCA += var
            ratio.append( ac/var )
        ac_mean.append( np.mean(ratio) )
        ac_std.append ( np.std(ratio) )         
            
        ac_TMRCA.append( cum_ac_TMRCA / cum_var_TMRCA )

    return ac_TMRCA, ac_mean, ac_std
    

def myfigures ( delta, rho, prefix, legend, colors):
    # rho is a list of MS, SCRM (pruned) and SCRM (full pruning) results
    # results is a list of autocorrelations, one for each delta
    #print legend
    l = [] 
    fig,(ax1)=pylab.subplots(1,1)
    for i in range ( len (rho) ):
        #y_err = [np.std(yi)*1.96/np.sqrt(len(yi)) for yi in rho[i]]
        tmp1 = ax1.plot( delta, rho[i] , color = colors[i])
        l.append ( tmp1 )
    pylab.xlim( [np.min(delta), np.max(delta)] )
    #pylab.title( prefix + " of " + `len(delta)` + " delta points" )
    pylab.xlabel(r'Distance between two sites $\delta$')
    pylab.ylabel(r'Autocorrelation $\rho$')
    pylab.legend ([ x[0] for x in l], legend, loc = 1)
    pylab.savefig( prefix+".pdf" )
    pylab.close()


def time_figure(accuracy, time, prefix, legend, colors):
    x = accuracy
    y = time
    markers = ["v", "o", "*", ">", "<", "s", "^", "+" , "D", "H"]
    #pylab.title("Time vs accuracy")
    pylab.ylabel("Time")
    pylab.xlabel("Accuracy")
    #pylab.xlabel("rho tmrca at delta = 10000")
    myl = []
    for i, xi in enumerate(x):
        myl.append(pylab.plot( x[i], np.log(y[i]), markers[i]))
    my_axes = pylab.gca()
    yticks = my_axes.get_yticks()
    ylabels = ["%.5g" % (np.exp(float(y))) for y in yticks]
    my_axes.set_yticklabels(ylabels)    
    pylab.legend( [ lx[0] for lx in myl ] ,  legend, loc=1, numpoints=1)
    pylab.savefig( prefix+"_timeVSacc.pdf")
    pylab.close()


def read_param_file ( experiment_name ):
    top_param = parameter()
    experiment_file = open( experiment_name, "r" )
    for line in experiment_file:
        if   line.split()[0] == "case:":          top_param.case       = line.split()[1]
        elif line.split()[0] == "nsam:":          top_param.nsam    = line.split()[1]
        elif line.split()[0] == "replicate:":     top_param.rep = int(line.split()[1])
        elif line.split()[0] == "seqlen:":        top_param.seqlen    = int(line.split()[1])
        elif line.split()[0] == "rho:":           top_param.rho = float(line.split()[1])
        elif line.split()[0] == "divergence:":    top_param.divergence    = int(line.split()[1])
        elif line.split()[0] == "job:":    top_param.jobs.append( line.split()[1] )
    experiment_file.close()
    
    top_param.printing()
    return top_param

def calculate_acurrcy(data_matrix, delta, obj_index =0 ): # obj_index is index for processed data, 0: tmrca, 1: TMRC, 2: clade
    accuracy = []
    ms_ld = data_matrix[0][obj_index]    
    print ms_ld
    for data_i in data_matrix:
        programs_ld = data_i[0]
        y = np.array([ ms_ld[i] - programs_ld[i] for i in range(len(ms_ld))] )
        accuracy.append(np.abs(simps(y, x = delta)))
    return accuracy

def calculate_acurrcy_array( data_array, delta): # obj_index is index for processed data, 0: tmrca, 1: TMRC, 2: clade
    accuracy = []
    ms_ld = data_array[0]
    for data_i in data_array:
        #programs_ld = data_i[0]
        y = np.array([ ms_ld[i] - data_i[i] for i in range(len(ms_ld))] )
        accuracy.append(np.abs(simps(y, x = delta)))
    return accuracy
    
if __name__ == "__main__":
    _use_param = read_param_file ( sys.argv[1] )
    _use_param.printing()
    
    for job in _use_param.jobs:
        print job
        data = extract_all_info (job, _use_param.rep)
        #print "len(data) = ", len(data)
        #print "len(data[0]) = ", len(data[0])
        #print data
        processed_data = process_data (data, _use_param.small_delta)
        print "sum over sum = ", processed_data[0]
        print "mean = ", processed_data[1]
        print "std = ", processed_data[2]
        print "mean and std of time " , process_time (data)
        #print "time = ", processed_data[3]
    
    #_legend = [ job[len(_use_param.case):-1] for job in _use_param.jobs ]
    #_colors = [ "orange", "purple", "green",  "red",  "blue", "black",  "yellow", "cyan", "magenta"]

    #for job_i, job in enumerate(_use_param.jobs):
        #print job_i
        #f = open ( job+"tmrcaRho", "w" )
        #f.write(`processed_data[job_i][0]`+"\n")
        #f.close()
        #f = open ( job+"time", "w" )
        #f.write(`processed_data[job_i][3]`+"\n")
        #f.close()
        ##print job

    
    #myfigures ( _use_param.small_delta, [ data_i[0] for data_i in processed_data ] , _use_param.case+"tmrca", _legend, _colors)
    #time_figure ( calculate_acurrcy(processed_data, _use_param.small_delta) , [ data_i[3] for data_i in processed_data ] , _use_param.case+"tmrca", _legend, _colors)
