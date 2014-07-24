#!/bin/bash
#$ -cwd
#$ -V
#$ -P bsg.prjb -q short.qb
#$ -e ErrFiles
#$ -o OutFiles
#$ -N fastsimcoal 
#$ -t 1-1000
#$ -j y


source parameters_preset


fsc_param_file=1Pop20sample.par

program=fastsimcoal

job=${case}${program}_

#for rep in $(seq 1 1 10)
    #do 
prefix=${job}${rep}
mkdir ${top_dir}"/"${prefix}
fileprefix=${top_dir}"/"${prefix}"/"${prefix}

infile=${prefix}.par
outfile=${prefix}"/"${prefix}_1_true_trees.trees    
cp ${fsc_param_file} ${infile}
echo ${fileprefix}
{ time -p ${program} -i ${infile} -n 1 -T --seed ${rep} > ${prefix}dummy ;} 2> ${fileprefix}timedummy.text
sed -e "/No/d" ${fileprefix}timedummy.text > ${fileprefix}time.text

grep ");" ${outfile} | sed -e "s/tree.*pos_/\\[/g" -e "s/ = \\[&U\\] /\\]/g" > ${fileprefix}

tree_file_name=${fileprefix}"Trees"
tree_change_name=${fileprefix}"change"
tree_freq_name=${fileprefix}"TreeFreq"
tmrca_raw_name=${fileprefix}"Tmrcaraw"
#tmrca_name=${fileprefix}"Tmrca"
first_coal_name=${fileprefix}"FirstCoal"

grep ';' ${fileprefix} | sed -e "s/\\[.*\\]//g" > ${tree_file_name}
grep ";" ${fileprefix} | sed -e "s/\\[//g" | sed -e "s/\\].*;//g" > ${tree_change_name}    
hybrid-Lambda -gt ${tree_file_name} -tmrca ${tmrca_raw_name}
#hybrid-Lambda -gt ${tree_file_name} -firstcoal ${first_coal_name}
./fastsimcoal_process.py ${fileprefix} 10000001

rm ${infile} ${outfile} ${fileprefix} ${tree_file_name} ${fileprefix}timedummy.text ${tree_change_name} ${tmrca_raw_name} ${prefix}dummy
rm -r ${prefix}

    
    #done
