#!/usr/bin/env python
import pylab as plt
from scipy.integrate import simps, trapz
import numpy as np

class job:
    def __init__(self, jobname):
        f = open ( jobname, "r" ) 
        self.ac = [float(x) for x in f.readline().strip('[').strip(']\n').split(',') ]
        self.ac_mean = [float(x) for x in f.readline().strip('[').strip(']\n').split(',') ]
        self.ac_std = [float(x) for x in f.readline().strip('[').strip(']\n').split(',') ]
        line = f.readline()
        self.time_mean = float(line.strip('[').strip(']\n').split(',')[0])
        self.time_std  = float(line.strip('[').strip(']\n').split(',')[1])
    
    def printing(self):
        print self.ac
        print self.ac_mean
        print self.ac_std
        print self.time_mean
        print self.time_std
        
    def diff_from_ms(self, ms_ac, delta):
        y = np.array([ np.abs(ms_ac[i] - self.ac[i]) for i in range(len(ms_ac))] )
        self.dev = np.abs(simps(y, x = delta))
        #print self.dev

delta = range( 0, int(2e4+1), 1000 )

prefix = "Constantpopsize"
suffix = "_Processed"
ms_case = [""]
fastsimcoal_case = [""]
#scrm_case = [ "window" + `x` for x in [0, 500, 1000, 3000, 5000, 7000, 10000, 30000, 50000, 70000, 100000] ]
#macs_case = [ "retain" + `x` for x in [0, 500, 1000, 3000, 5000, 7000, 10000, 30000, 50000, 70000, 100000] ]
scrm_case = [ "window" + `x` for x in [0,  5000, 7000, 10000, 30000, 50000, 100000] ]
macs_case = [ "retain" + `x` for x in [0,  5000, 7000, 10000, 30000, 50000, 100000] ]
program = ["ms", "fastsimcoal", "scrm", "macs"]
case = [ms_case, fastsimcoal_case, scrm_case, macs_case]
joblist = []

ms = job("Constantpopsizems_Processed")
#ms.printing()

fig1 = plt.figure(figsize=(9, 6), dpi=80)
ax1 = fig1.add_subplot(111)

fig2 = plt.figure(figsize=(12, 8), dpi=80)
ax2 = fig2.add_subplot(111)

fig3 = plt.figure(figsize=(9, 6), dpi=80)
ax3 = fig3.add_subplot(111)

#ax1.plot ( delta, ms.ac , linewidth=3.0, color = "black")
#ax1.errorbar ( delta, ms.ac, yerr = [ x/(1000**0.5) *1.96 for x in ms.ac_std ] )

linestyles = ['-', '-.', '--', ':']
colors = [ "blue", "red", "green",  "purple", "black",  "yellow", "cyan", "magenta", "orange"]
markers = ["v", "o", "*", ">", "<", "s", "^", "+" , "D", "H", "d","x"]
legendlist1 = []
l1 = []

legendlist2 = []
l2 = []

legendlist3 = []
l3 = []

for i, program_i in enumerate ( program ):
    color_j = 0
    program_dev = []
    program_time = []
    program_time_err = []
    for j, case_j in enumerate ( case[i] ):
        current_job2 = job( prefix + program_i + case_j + suffix )
        current_job2.diff_from_ms (ms.ac, delta) 
        program_dev.append ( current_job2.dev)
        program_time.append (current_job2.time_mean)
        program_time_err.append( current_job2.time_std  ) 
        #current_dot = ax2.plot ( current_job2.dev, np.log(current_job2.time_mean), markers[j], color = colors[i])
        current_dot = ax2.plot ( current_job2.dev, current_job2.time_mean, markers[j], color = colors[i])
        l2.append(current_dot)
        legendlist2.append(program_i + case_j)
        if j % 2 == 0:
            current_job = job( prefix + program_i + case_j + suffix )
            legendlist1.append( program_i + case_j)
            current_line = ax1.plot ( delta, current_job.ac, linestyles[i], color = colors[color_j] )
#           ax1.errorbar ( delta, current_job.ac, yerr = [ x/(1000**0.5) *1.96 for x in current_job.ac_std ],
#                        fmt='.', color = colors[color_j] )
            l1.append(current_line)
            relative_ac = [ np.abs(ms.ac[ac_i] - current_job.ac[ac_i]) for ac_i in range(len(ms.ac))]
            current_line3 = ax3.plot ( delta, relative_ac, linestyles[i], color = colors[color_j] )
 #           ax3.errorbar ( delta, relative_ac, yerr = [ x/(1000**0.5) *1.96 for x in current_job.ac_std ],
 #                         fmt='.', color = colors[color_j] )
            l3.append(current_line3)
            color_j += 1
            
    #ax2.errorbar ( program_dev, np.log(program_time), yerr = program_time_err, color = colors[i])
    ax2.plot ( program_dev, program_time, color = colors[i])
    #ax2.errorbar ( program_dev, program_time, yerr = program_time_err, color = colors[i])
    
ms_line = ax1.plot ( delta, ms.ac,  linewidth=2.0, color = "black")    
#ax1.errorbar ( delta, ms.ac, yerr = [ x/(1000**0.5) *1.96 for x in current_job.ac_std ],
                          #fmt='.', color = colors[color_j] )
l1[0] = ms_line                          

relative_ac = [ float(0) for ac_i in range(len(ms.ac))]
ms_line3 = ax3.plot ( delta, relative_ac,  linewidth=2.0, color = "black")    
#ax3.errorbar ( delta, relative_ac, yerr = [ x/(1000**0.5) *1.96 for x in ms.ac_std ],
                          #fmt='.', color = "black" )
l3[0] = ms_line3
                              
ax1.legend ([ x[0] for x in l1], legendlist1, loc = 1)        
ax1.axis([0,20000, 0, 1]) 
ax1.set_xlabel(r'Distance between two sites $\delta$')
ax1.set_ylabel(r'Autocorrelation $\rho$')
fig1.savefig("TMRCArhoLD.pdf")

ax3.legend ([ x[0] for x in l1], legendlist1, loc = 1)        
ax3.axis([0,30000, -.01, 0.06]) 
ax3.set_xlabel(r'Distance between two sites $\delta$')
ax3.set_ylabel(r'Error in Autocorrelation $\rho$')
fig3.savefig("RelativeTMRCArhoLD.pdf")


ax2.set_xlim ([-10, 1300])
ax2.set_ylim ([0.5, 35])
#ax2.axis([-10, 1300, -2.5, np.log(ms.time_mean)*1.1] )
#yticks = ax2.get_yticks()
#print yticks
#ylabels = ["%.5g" % (np.exp(float(y))) for y in yticks]
#ax2.set_yticklabels(ylabels) 
ax2.set_yscale('log')
ax2.legend ([ x[0] for x in l2], legendlist2, loc = 1, numpoints=1)        
ax2.set_ylabel("Time (sec)")
ax2.set_xlabel("Deviation")
fig2.savefig("time_vs_dev.pdf")
