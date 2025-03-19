#! /bin/bash 
#PBS -l nodes=1:ppn=1:centos7,cput=24:00:00,walltime=48:00:00 
#PBS -N process_data
#PBS -d /export/biostuds/2266643a/BIOL5177_rnaseq/biol5177_work/cw_1
#PBS -m abe 
#PBS -M 2266643a@student.gla.ac.uk 
#PBS -q bioinf-stud

# RESOURCES
adapter=illumina_adapter.fa # path to adapter
hs2index="/export/projects/polyomics/Genome/Mus_musculus/mm10/Hisat2Index/chr2"         # path to reference genome hisat2 indexes 
gtf=/export/projects/polyomics/Genome/Mus_musculus/mm10/annotations/chr2.gtf   	         # path to GTF file 
data=./data             # path to local data directory that you have created 

# <Make subdirs
hisat_dir=./hisat_results
stringtie_dir=./stringtie_results
mkdir -p $hisat_dir
mkdir -p $stringtie_dir

gtflist='list.gtf.txt' # filename for a final GTF list 

rm -f ${gtflist}            # remove if exists 

# Running in loop to carry out work
for sample in s1 s1t s2 s2t s3 s4 s5 s6 s7 s8 s9 s10 s11 s12

do
    sample="${sample}.c2.25K"
    fastq="${data}/${sample}.fq"	# path to raw fastq file 
    trim1="$data/${sample}.t1.fq"		# path to adapter-trimmed fastq file 
    trim2="$data/${sample}.t2.fq"		# path to quality-trimmed fastq file 
    sam="${hisat_dir}/${sample}.sam"  # path to hisat2-generated SAM file
    bam="${hisat_dir}/${sample}.bam"  # path to samtools-converted BAM file 
    sorted_bam="${hisat_dir}/${sample}.sort.bam"	 # path to samtools-sorted BAM file    
    
    
    scythe -o $trim1 -a illumina_adapter.fa -q sanger ${fastq}
    sickle se -f $trim1 -t sanger -o $trim2  -q 10 -l 51
    hisat2 -p 4 --phred33 -x ${hs2index} -U ${trim2} -S ${sam} # command to execute hisat2 
    samtools view -b -o ${bam} ${sam}		# command to execute samtools view 
    samtools sort -o ${sorted_bam} ${bam}	# command to execute samtools sort 

    rm ${sam} ${bam}		# removing unnecessary SAM/BAM files 
    rm ${trim1} ${trim2}	# removing unnecessary trimmed FASTQ files 

    str_smp_dir="${stringtie_dir}/${sample}"		# path to sample-specific subdirectory for stringtie results 

    mkdir -p $str_smp_dir	# make the above directory 
    sample_tr_gtf="${str_smp_dir}/${sample}_transcripts.gtf"
    stringtie -p 4 -t -e -B -G ${gtf} -o ${sample_tr_gtf} ${sorted_bam}	# command to execute stringtie 

    gtfline="${sample} $sample_tr_gtf" # line containing sample and path to string tie-generated GTF file 
    echo ${gtfline} >> ${gtflist}	# adding the above line to a file 

done 

prepDE.py -i ${gtflist}
