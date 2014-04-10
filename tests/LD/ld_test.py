#!/usr/bin/env python

import pylab as pl
import numpy as np
import os
import sys
import pylab

__mydebug__       = False
__fix_ms_seed__   = False
__fix_scrm_seed__ = False

# exact_windows_length [0, 10e2, 10e3, 10e4, -1]

class parameter:
    def __init__( self, nsam = 6, replicate = 100, seqlen = 1e7, rho = 1e-8 , exact_window_length = 0 , divergence = False): 
        """
        Define simulation parameters
        """
        self.nsam      = nsam
        self.rep       = replicate
        self.seqlen    = seqlen # 10e7
        #self.theta     = 7 * 10e-10 * seqlen * 4 * 10000
        self.rho       = rho * (self.seqlen-1) * 4 * 10000
        self.exact_window_length = exact_window_length
        self.divergence = divergence
        
        big_delta_max = 1e5
        small_delta_max = 1e4
        self.big_delta = range( 0, int(big_delta_max+1), 10000 )
        self.small_delta = range( 0, int(small_delta_max+1), 1000 )

    
    def printing ( self ):
        """
        Check simulation parameters
        """
        print "sample size:       ", self.nsam
        print "replicate:         ",   self.replicate
        print "Sequence length:   ", self.seqlen
        print "recombination rate:", self.rho
        
    
    def define_command ( self, scrm = False ):        
        cmd = `self.nsam` + " 1 " + " -T " + " -r " + `self.rho` + " " + `int(self.seqlen)`
        if self.divergence: cmd += " -I 2 " + `self.nsam/2` + " " + `self.nsam/2` + " -ej 1 2 1 "
        if scrm :  cmd += " -l " + `int(self.exact_window_length)`
        if __fix_ms_seed__ :   cmd += " -seed 2 2 2 "
        if __fix_scrm_seed__ :   cmd += " -seed 2 "
        if __mydebug__: print cmd
        #print cmd
        return cmd
        
        
        
def run_program ( param, scrm = False, prefix = "msout" ):
    program = "{ time -p "
    program += "ms" if not scrm else "scrm"
    program += " " + param.define_command( scrm ) + " > " + prefix + " ;} 2> "+prefix+"time.text"
    #print program        
    os.system ( program )    
    return prefix
    
    
def process_ms_scrm_output ( prefix ):
    tree_file_name  = prefix + "Trees"
    tree_freq_name  = prefix + "TreeFreq"
    tmrca_name      = prefix + "Tmrca"
    first_coal_name = prefix + "FirstCoal"
    
    grep_tree       = "grep \';\' " + prefix + " | sed -e 's/\\[.*\\]//g' > " + tree_file_name  
    grep_freq       = "grep \';\' " + prefix + " | sed -e 's/\\[//g' | sed -e 's/\\].*;//g' > " + tree_freq_name
    grep_tmrca      = "hybrid-Lambda -gt " + tree_file_name + " -log -tmrca " + tmrca_name
    grep_first_coal = "hybrid-Lambda -gt " + tree_file_name + " -log -firstcoal " + first_coal_name
    
    os.system( grep_tree       )
    os.system( grep_freq       )
    os.system( grep_tmrca      )
    os.system( grep_first_coal )
    
    tree_freq = [ float(x) for x in open( tree_freq_name, "r" ) ]
    
    tmrca = [ float(x) for x in open( tmrca_name, "r" )]
    
    first_coal_file = open( first_coal_name, "r" )
    first_coal_time = []
    clade = []
    for line in first_coal_file:
        first_coal_time.append( float(line.split()[0]) )        
        clade.append( line.split()[1] )
    first_coal_file.close()
    
    timeFile_name = prefix+"time.text"
    runtime = float(open( timeFile_name, "r").readlines()[1].strip("user").strip('\n'))
    #print runtime

    return  tree_freq, tmrca, first_coal_time, clade, runtime


def cal_ac_TMRC_star (tree_freq, tmrca, avg_tmrca, delta):

    seqlen = sum(tree_freq)

    cumfreq = [0]
    shifted_cumfreq = [] 
    for x in tree_freq : 
        cumfreq.append(cumfreq[-1] + x)
        shifted_cumfreq.append(cumfreq[-1] - delta)
    cumfreq.pop(0)

    n = int(seqlen - delta)
    ac = 0
    var = 0
    
    if __mydebug__:
        print tree_freq
        print cumfreq
        print shifted_cumfreq        
        print "seqence length = ", seqlen
        print "i should iterate until ", n, "trees"
    
    T1_index = 0
    T2_index = 0 # use while loop to determine the initial T2_index
    while delta > cumfreq[T2_index]: # need to check again ...
        T2_index += 1
    about_to = False
    for i in range( n ): 
        if __mydebug__ & ((not i < cumfreq[T1_index]) | (not i < shifted_cumfreq[T2_index])):
            print "before", i, T1_index,  T2_index
            print "-------------------------------"
            about_to = True
        if i >= cumfreq[T1_index]: T1_index += 1
        term1 = tmrca[T1_index] - avg_tmrca
        if i >= shifted_cumfreq[T2_index]: T2_index += 1
        term2 = tmrca[T2_index] - avg_tmrca
        if __mydebug__ & about_to:
            print "after ", i, T1_index,  T2_index
            about_to = False
        #print i, T1_index, term1, T2_index, term2
        ac += term1 * term2
        var += term1 * term1
    ac /= float(n)
    var /= float(n)
    return ac, var


def cal_ac_TMRC_star_2 (tree_freq, tmrca, avg_tmrca, delta):

    seqlen = sum(tree_freq)

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


def cal_ac_clade_2 (tree_freq, clade, delta ):
    seqlen = sum(tree_freq)
    cumfreq = [0]
    for x in tree_freq : 
        cumfreq.append(cumfreq[-1] + x)
    cumfreq.pop(0)
    n = int(seqlen - delta)
    ac = 0.0
    length = 0.0
    T1_index = 0
    T2_index = 0 # use while loop to determine the initial T2_index
    while delta > cumfreq[T2_index]: # need to check again ...
        T2_index += 1
    i = 0
    while i < n:
        distance = min( cumfreq[T1_index] - i, 
                        cumfreq[T2_index] - delta - i )
        distance = min( distance, n - i )
        term1 = clade[T1_index]  
        term2 = clade[T2_index]  
        ac += distance if term1 == term2 else 0
        length += distance
        i += distance            
        if i >= cumfreq[T1_index]: T1_index += 1
        if i + delta >= cumfreq[T2_index]: T2_index += 1

    ac /= float(n)
    length /= float(n)   # this should be 1
    return ac, length

    
    
def cal_ac_clade (tree_freq, clade, delta ):
    seqlen = sum(tree_freq)
    cumfreq = [0]
    shifted_cumfreq = [] 
    for x in tree_freq : 
        cumfreq.append(cumfreq[-1] + x)
        shifted_cumfreq.append(cumfreq[-1] - delta)
    cumfreq.pop(0)
    #print cumfreq, shifted_cumfreq
    n = int(seqlen - delta)
    ac = 0.0
    length = 0.0
    T1_index = 0
    T2_index = 0 # use while loop to determine the initial T2_index
    while delta > cumfreq[T2_index]: # need to check again ...
        T2_index += 1
    for i in range( n ):        
        T1_index += 0 if i < cumfreq[T1_index] else 1
        term1 = clade[T1_index] 
        T2_index += 0 if i < shifted_cumfreq[T2_index] else 1
        term2 = clade[T2_index] 
        #print i, T1_index, term1, T2_index, term2
        ac += 1 if term1 == term2 else 0
        length += 1
    ac /= float(n)
    length /= float(n)   # this should be 1
    return ac, length
    
    
def n_rep ( param, scrm, prefix ):

    reps = param.rep

    # first collect all the data
    data = []
    for rep in range(reps):
        if rep%(reps/10) == 0: print "replicate:",rep
        out = run_program ( param, scrm, prefix)
        # (tree_freq, tmrca, first_coal_time, clade, runtime)
        data.append( process_ms_scrm_output ( out ) )
        
    # compute the average Tmrca and Tmrc
    tot_tmrca = 0
    tot_tmrc = 0
    tot_time = 0
    tot_runtime = 0
    for d in data:
        tot_runtime += d[4]
        for i, duration_i in enumerate( d[0] ):
            tmrca_i = d[1][i]
            tmrc_i = d[2][i]
            tot_tmrca += tmrca_i
            tot_tmrc += tmrc_i
            tot_time += duration_i
    avg_tmrca = tot_tmrca / tot_time
    avg_tmrc = tot_tmrc / tot_time

    ac_TMRC = []    
    ac_clade = []
    for delta_i in param.big_delta:
        cum_ac_TMRC = 0
        cum_var_TMRC = 0
        cum_ac_clade = 0
        cum_length_clade = 0
        for d in data:
            ac, var = cal_ac_TMRC_star_2( d[0], d[2], avg_tmrc, delta_i )
            cum_ac_TMRC += ac
            cum_var_TMRC += var
            ac, length = cal_ac_clade_2( d[0], d[3], delta_i )
            cum_ac_clade += ac
            cum_length_clade += length
        ac_TMRC.append( cum_ac_TMRC / cum_var_TMRC )
        ac_clade.append( cum_ac_clade / cum_length_clade )
        
    ac_TMRCA  = []    
    for delta_i in param.small_delta:        
        cum_ac_TMRCA = 0
        cum_var_TMRCA = 0
        for d in data:
            ac, var = cal_ac_TMRC_star_2( d[0], d[1], avg_tmrc, delta_i )
            cum_ac_TMRCA += ac
            cum_var_TMRCA += var
        ac_TMRCA.append( cum_ac_TMRCA / cum_var_TMRCA )

    return ac_TMRCA, ac_TMRC, ac_clade, tot_runtime
    

#def n_rep ( param, scrm , prefix):
    #ac = []
    
    #for i in range( param.rep ):
        #if i%10 == 0: print i
        #ac.append ( one_single_rep ( param, scrm, prefix ) )
    ## ac[0] = ac_TMRCA, which has length of 10
    
    #Etmrca = [] # two dimension [i][j] i is ith replication, j is the Etmrca associated with delta[j]
    #Efirst_coal = []
    #Eclade = []
    
    #time = [ x[3] for x in ac ]
    
    #for i in range( len(param.big_delta) ):
        #Etmrca.append( [ x[0][i] for x in ac ] ) # x is one single repicate, x[0] is the ac_TMRCA, x[0][i] is the ith Etmrca associated with delta[i]
        #Efirst_coal.append( [ x[1][i] for x in ac ] )
        #Eclade.append ( [ x[2][i] for x in ac ] )
        
    #return Etmrca, Efirst_coal, Eclade, time


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
    pylab.xlabel(r'Distance between two sites $\delta$')
    pylab.ylabel(r'Autocorrelation $\rho$')
    pylab.legend ([ x[0] for x in l], legend, loc = 1)
    pylab.savefig( prefix+".pdf" )
    pylab.close()


if __name__ == "__main__":

    _use_param = parameter(nsam = 6, replicate = 1000, seqlen = 1e7, rho = 1e-8)         
    _msac = n_rep ( _use_param, scrm = False, prefix = "msout")        
    _scrmac = n_rep ( _use_param, scrm = True, prefix = "scrmout") 
    
    _use_param.exact_window_length = 1e3
    _scrmace3 = n_rep ( _use_param, scrm = True, prefix = "scrme3out")        

    _use_param.exact_window_length = 1e4
    _scrmace4 = n_rep ( _use_param, scrm = True, prefix = "scrme4out")      

    _use_param.exact_window_length = 5e4
    _scrmac5e4 = n_rep ( _use_param, scrm = True, prefix = "scrm5e4out")      

    _use_param.exact_window_length = 1e5
    _scrmace5 = n_rep ( _use_param, scrm = True, prefix = "scrme5out")        
    
    _legend = ["ms", "exact window = 10000", "exact window = 1000", "exact window = 0"]
    _colors = ["red", "green", "blue", "black"]
    
    # extract TMRCA from results and plot
    _rho = [ _msac[0], _scrmace5[0], _scrmace3[0], _scrmac[0] ]
    myfigures ( _use_param.small_delta, _rho, "tmrca", _legend, _colors)
    
    # extract TMRC
    _rho = [ _msac[1], _scrmace5[1], _scrmace3[1], _scrmac[1] ]
    myfigures ( _use_param.big_delta, _rho, "tmrc", _legend, _colors)

    # extract clade
    _rho = [ _msac[2], _scrmace5[2], _scrmace3[2], _scrmac[2] ]
    myfigures ( _use_param.big_delta, _rho, "clade", _legend, _colors)
    
    
    ########################## Divergence
    _use_paramD = parameter(nsam = 6, replicate = 1000, seqlen = 1e7, rho = 1e-8, divergence = True)         
    _msacD = n_rep ( _use_paramD, scrm = False, prefix = "Divergencemsout")        
    _scrmacD = n_rep ( _use_paramD, scrm = True, prefix = "Divergencescrmout") 
    
    _use_param.exact_window_length = 1e3
    _scrmace3D = n_rep ( _use_paramD, scrm = True, prefix = "Divergencescrme3out")        

    _use_param.exact_window_length = 1e4
    _scrmace4D = n_rep ( _use_paramD, scrm = True, prefix = "Divergencescrme4out")        

    _use_param.exact_window_length = 5e4
    _scrmac5e4D = n_rep ( _use_paramD, scrm = True, prefix = "Divergencescrm5e4out")        

    _use_param.exact_window_length = 1e5
    _scrmace5D = n_rep ( _use_paramD, scrm = True, prefix = "Divergencescrme5out")        
    
    
    # extract TMRCA from results and plot
    _rho = [ _msacD[0], _scrmace5D[0], _scrmace3D[0], _scrmacD[0] ]
    myfigures ( _use_param.small_delta, _rho, "Divergencetmrca", _legend, _colors)
    
    # extract TMRC
    _rho = [ _msacD[1], _scrmace5D[1], _scrmace3D[1], _scrmacD[1] ]
    myfigures ( _use_param.big_delta, _rho, "Divergencetmrc", _legend, _colors)

    # extract clade
    _rho = [ _msacD[2], _scrmace5D[2], _scrmace3D[2], _scrmacD[2] ]
    myfigures ( _use_param.big_delta, _rho, "Divergenceclade", _legend, _colors)
    #except:
        ##print "Usage: %s  <seqlen>  <position_file_name>  <psmc_input_file_prefix>" % sys.argv[0]
        #sys.exit(1)
