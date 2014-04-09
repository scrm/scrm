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
    def __init__( self, nsam = 5, replicate = 1, seqlen = 1e7, rho = 1e-8 , exact_window_length = 0 ): 
        """
        Define simulation parameters
        """
        self.nsam      = nsam
        self.replicate = replicate
        self.seqlen    = seqlen # 10e7
        #self.theta     = 7 * 10e-10 * seqlen * 4 * 10000
        self.rho       = rho * (self.seqlen-1) * 4 * 10000
        self.exact_window_length = exact_window_length
        self.big_delta = [100, 1000, 2000]
        #self.small_delta = 
    
    def printing ( self ):
        """
        Check simulation parameters
        """
        print "sample size:       ", self.nsam
        print "replicate:         ",   self.replicate
        print "Sequence length:   ", self.seqlen
        print "recombination rate:", self.rho
        
    
    def define_command ( self, scrm = False ):        
        cmd = `self.nsam` + " " + `self.replicate` + " -T " + " -r " + `self.rho` + " " + `self.seqlen`
        if scrm :  cmd += " -l " + `self.exact_window_length`            
        if __fix_ms_seed__ :   cmd += " -seed 2 2 2 "
        if __fix_scrm_seed__ :   cmd += " -seed 2 "
        if __mydebug__: print cmd
        #print cmd
        return cmd
        
        
def run_program ( param, scrm = False, prefix = "msout" ):
    program = "{ time -p "
    program += "ms" if not scrm else "scrm"
    program += " " + param.define_command( scrm ) + " > " + prefix + " ;} 2> time.text"
    print program        
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
    
    runtime = float(open("time.text", "r").readlines()[1].strip("user").strip('\n'))
    #print runtime

    return  tree_freq, tmrca, first_coal_time, clade, runtime



        
    
#def extract_info ( outfiles ):
    #tree_freq_file = open( outfiles[0], "r" )
    #tree_freq = [ float(x) for x in tree_freq_file]
    #tree_freq_file.close()
    ##total_num_tree = 0
    ##for freq in tree_freq : total_num_tree += freq
    ##normalized_freq = [ x / float(total_num_tree) for x in tree_freq ]
    ##print normalized_freq
    
    #tmrca_file = open( outfiles[1], "r" )
    #tmrca = [ float(x) for x in tmrca_file]
    #tmrca_file.close()
    
    #first_coal_file = open( outfiles[2], "r" )
    #first_coal_time = []
    #clade = []
    #for line in first_coal_file:
        #first_coal_time.append( float(line.split()[0]) )        
        #clade.append( line.split()[1] )
    #first_coal_file.close()
    
    #return  tree_freq, tmrca, first_coal_time, clade


def cal_ac_TMRC_star (tree_freq, tmrca, param, delta):
    seqlen = param.seqlen    
    T_bar = 0
    for i in range( len( tree_freq ) ):
        T_bar = tree_freq[i] * tmrca[i]
    T_bar /= float(seqlen)
    
    #print T_bar
    
    cumfreq = [0]
    shifted_cumfreq = [] 
    for x in tree_freq : 
        cumfreq.append(cumfreq[-1] + x)
        shifted_cumfreq.append(cumfreq[-1] - delta)
    cumfreq.pop(0)

    n = int(seqlen - delta)
    ac = 0
    
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
        T1_index += 0 if i < cumfreq[T1_index] else 1 
        term1 = tmrca[T1_index] - T_bar
        T2_index += 0 if i < shifted_cumfreq[T2_index] else 1 
        term2 = tmrca[T2_index] - T_bar
        if __mydebug__ & about_to:
            print "after ", i, T1_index,  T2_index
            about_to = False
        #print i, T1_index, term1, T2_index, term2
        ac += term1 * term2
    ac /= float(n)
    return ac
    
    
def cal_ac_clade (tree_freq, clade, param, delta ):
    seqlen = param.seqlen       
    cumfreq = [0]
    shifted_cumfreq = [] 
    for x in tree_freq : 
        cumfreq.append(cumfreq[-1] + x)
        shifted_cumfreq.append(cumfreq[-1] - delta)
    cumfreq.pop(0)
    #print cumfreq, shifted_cumfreq
    n = int(seqlen - delta)
    ac = 0.0
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
    ac /= float(n)
    return ac
    
    
def one_single_rep ( param, scrm, prefix ):
    out = run_program ( param, scrm, prefix)
    (tree_freq, tmrca, first_coal_time, clade, runtime) = process_ms_scrm_output ( out )
    #info = extract_info ( outfiles )
        
    ac_TMRCA = []
    ac_TMRC  = []
    ac_clade = []
    for delta_i in param.big_delta:
        ac_TMRCA.append( cal_ac_TMRC_star ( tree_freq, tmrca, param, delta_i ) )        
        ac_clade.append( cal_ac_clade ( tree_freq, clade, param, delta_i ) )
        
        #ac_TMRC.append ( cal_ac_TMRC_star ( tree_freq, first_coal_time, param, delta_i ) )
    return ac_TMRCA, ac_clade
    
    #return  

def n_rep ( param, delta, scrm , prefix, rep = 10000):
    ac = []
    for i in range(rep):
        ac.append ( one_single_rep ( param, delta, scrm, prefix ) )
    Etmrca = [ x[0] for x in ac ]
    Efirst_coal = [ x[1] for x in ac ]
    Eclade = [ x[2] for x in ac ]
    return Etmrca, Efirst_coal, Eclade

if __name__ == "__main__":
    try:
        use_param = parameter()
        #msac = n_rep (use_param, delta = 10000, scrm = False, prefix = "msout")
        #scrmac = n_rep (use_param, delta = 10000, scrm = True, prefix = "scrmout")

    except:
        #print "Usage: %s  <seqlen>  <position_file_name>  <psmc_input_file_prefix>" % sys.argv[0]
        sys.exit(1)
