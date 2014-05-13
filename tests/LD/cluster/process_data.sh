#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e ErrFiles
#$ -o OutFiles
#$ -N process_data.sh
#$ -t 1-3
#$ -j y

case="Constantpopsize"

joblist=("ms"  "scrm" "macs")

rep=0 # starts from zero
echo ${joblist[${rep}]}



#case: Constantpopsize
#nsam: 6
#replicate: 1000
#seqlen: 10000001
#rho: 4000
#job: Constantpopsizems_
#job: Constantpopsizescrmwindow100000_
#job: Constantpopsizescrmwindow70000_
#job: Constantpopsizescrmwindow50000_
#job: Constantpopsizescrmwindow30000_
#job: Constantpopsizescrmwindow10000_
#job: Constantpopsizescrmwindow7000_
#job: Constantpopsizescrmwindow5000_
#job: Constantpopsizescrmwindow3000_
#job: Constantpopsizescrmwindow1000_
#job: Constantpopsizescrmwindow500_
#job: Constantpopsizescrmwindow0_
#job: Constantpopsizefastsimcoal_
#job: Constantpopsizemacs_
#job: Constantpopsizemacsretain1000_
#job: Constantpopsizemacsretain10000_
#job: Constantpopsizemacsretain100000_
