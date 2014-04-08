#!/usr/bin/env python

import pylab as pl
import numpy as np
import os
import sys

# exact_windows_length [0, 10e2, 10e3, 10e4, -1]

class parameter:
    def __init__( self, nsam = 20, replicate = 1, seqlen = 1000000, rho = 8 * 10e-10 , exact_window_length = 0 ): 
        self.nsam      = nsam
        self.replicate = replicate
        self.seqlen    = seqlen # 10e7
        #self.theta     = 7 * 10e-10 * seqlen * 4 * 10000
        self.rho       = rho * (self.seqlen-1) * 4 * 10000
        self.exact_window_length = exact_window_length
        
    def define_command ( self, scrm = False ):        
        cmd = `self.nsam` + " " + `self.replicate` + " -T " + " -r " + `self.rho` + " " + `self.seqlen`
        if scrm :
            cmd += " -l " + `self.exact_window_length`
        return cmd


def process_ms_scrm_output ( prefix ):
    tree_file  = prefix + "Trees"
    tree_freq  = prefix + "TreeFreq"
    tmrca      = prefix + "Tmrca"
    first_coal = prefix + "FirstCoal"
    
    grep_tree       = "grep \';\' " + prefix + " | sed -e 's/\\[.*\\]//g' > " + tree_file  
    grep_freq       = "grep \';\' " + prefix + " | sed -e 's/\\[//g' | sed -e 's/\\].*;//g' > " + tree_freq
    grep_tmrca      = "hybrid-Lambda -gt " + tree_file + " -log -tmrca " + tmrca
    grep_first_coal = "hybrid-Lambda -gt " + tree_file + " -log -firstcoal " + first_coal

    os.system( grep_tree       )
    os.system( grep_freq       )
    os.system( grep_tmrca      )
    os.system( grep_first_coal )
    return tree_freq, tmrca, first_coal


def run_ms ( param, prefix = "msout" ):
    ms_command = "ms " + param.define_command() + " > " + prefix
    print ms_command        
    os.system ( ms_command )
    return prefix
    
    
def run_scrm ( param, prefix = "scrmout", scrm = False ):
    scrm_command = "scrm " + param.define_command( scrm ) + " > " + prefix
    print scrm_command        
    os.system ( scrm_command )
    return prefix
    

def extract_info ( outfiles ):

    tree_freq_file = open( outfiles[0], "r" )
    tree_freq = [ float(x) for x in tree_freq_file]
    tree_freq_file.close()
    #total_num_tree = 0
    #for freq in tree_freq : total_num_tree += freq
    #normalized_freq = [ x / float(total_num_tree) for x in tree_freq ]
    #print normalized_freq
    
    tmrca_file = open( outfiles[1], "r" )
    tmrca = [ float(x) for x in tmrca_file]
    tmrca_file.close()
    
    first_coal_file = open( outfiles[2], "r" )
    first_coal_time = []
    clade = []
    for line in first_coal_file:
        first_coal_time.append( float(line.split()[0]) )        
        clade.append( line.split()[1] )
    first_coal_file.close()
    
    return  tree_freq, tmrca, first_coal_time, clade


def ac_TMRC_star (tree_freq, tmrca, param, delta = 100000):
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
    #print cumfreq, shifted_cumfreq
    n = seqlen - delta
    ac = 0
    #print "seqence length = ", seqlen
    #print "i should iterate until ", n, "trees"
    T1_index = 0
    T2_index = 0
    for i in range( n ):        
        T1_index += 0 if i < cumfreq[T1_index] else 1
        term1 = tmrca[T1_index] - T_bar
        T2_index += 0 if i < shifted_cumfreq[T2_index] else 1
        term2 = tmrca[T2_index] - T_bar
        #print i, T1_index, term1, T2_index, term2
        ac += term1 * term2
    #ac /= float(n)
    return ac
    
    
def ac_clade (tree_freq, clade, param, delta = 100000):
    seqlen = param.seqlen    
   
    cumfreq = [0]
    shifted_cumfreq = [] 
    for x in tree_freq : 
        cumfreq.append(cumfreq[-1] + x)
        shifted_cumfreq.append(cumfreq[-1] - delta)
    cumfreq.pop(0)
    #print cumfreq, shifted_cumfreq
    n = seqlen - delta
    ac = 0
    #print "seqence length = ", seqlen
    #print "i should iterate until ", n, "trees"
    T1_index = 0
    T2_index = 0
    for i in range( n ):        
        T1_index += 0 if i < cumfreq[T1_index] else 1
        term1 = clade[T1_index] 
        T2_index += 0 if i < shifted_cumfreq[T2_index] else 1
        term2 = clade[T2_index] 
        #print i, T1_index, term1, T2_index, term2
        ac += 1 if term1 == term2 else 0
    #ac /= float(n)
    return ac
    
    
def one_single_rep ( param, delta = 10000, scrm = False ):
    out = run_ms ( param )
    outfiles = process_ms_scrm_output ( out )
    info = extract_info ( outfiles )
    
    #print ac_TMRC_star (info[0], info[1], use_param, delta = 10000)
    #print ac_TMRC_star (info[0], info[2], use_param, delta = 10000)
    #print ac_clade     (info[0], info[3], use_param, delta = 10000)
    return  ac_TMRC_star (info[0], info[1], param, delta ), \
            ac_TMRC_star (info[0], info[2], param, delta ), \
            ac_clade     (info[0], info[3], param, delta )

if __name__ == "__main__":
    try:
        use_param = parameter()
        ac = one_single_rep (use_param)
        #scrmout = run_scrm ( use_param, scrm = True)
        #scrmoutfiles = process_ms_scrm_output ( scrmout )

    except:
        #print "Usage: %s  <seqlen>  <position_file_name>  <psmc_input_file_prefix>" % sys.argv[0]
        sys.exit(1)
