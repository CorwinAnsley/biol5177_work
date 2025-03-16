#! /bin/bash 
#PBS -l nodes=1:ppn=4,cput=24:00:00,walltime=48:00:00 
#PBS -N link_data 
#PBS -d /export/biostuds/2266643a/BIOL5177_rnaseq/biol5177_work/cw_1
#PBS -m abe 
#PBS -M 2266643a@student.gla.ac.uk 
#PBS -q bioinf-stud

fqpath=/export/home/gmh5n/assignment
for sample in s1 s1t s2 s2t s3 s4 s5 s6 s7 s8 s9 s10 s11 s12

do
    raw_file=$fqpath/${sample}.fq  # path to the raw fast file 
    link_name=/export/biostuds/2266643a/BIOL5177_rnaseq/biol5177_work/cw_1/data/${sample}.fq     # link-name to appear in your directory 

    ln -s $raw_file $link_name             # instruction to make a soft link 
    test_file=/export/biostuds/2266643a/BIOL5177_rnaseq/biol5177_work/cw_1/data/${sample}.25K.fq       # test fast file containing 25K reads 
    head -n 100000  $link_name > $test_file       # instruction to generate test file 
done 
