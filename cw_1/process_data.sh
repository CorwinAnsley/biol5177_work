#! /bin/bash 
#PBS -l nodes=1:ppn=1:centos7,cput=24:00:00,walltime=48:00:00 
#PBS -N process_data
#PBS -d /export/biostuds/2266643a/BIOL5177_rnaseq/biol5177_work/cw_1
#PBS -m abe 
#PBS -M 2266643a@student.gla.ac.uk 
#PBS -q bioinf-stud

# RESOURCES - paths to various files need to run the cript
adapter=illumina_adapter.fa 
hs2index="/export/projects/polyomics/Genome/Mus_musculus/mm10/Hisat2Index/chr2"         
gtf=/export/projects/polyomics/Genome/Mus_musculus/mm10/annotations/chr2.gtf   	    
data=./data             

# Make the necessary subdirs
hisat_dir=./hisat_results
stringtie_dir=./stringtie_results
mkdir -p $hisat_dir
mkdir -p $stringtie_dir

# Define filename for gtf list
gtflist='list.gtf.txt' 

# Remove gtflist if exists already
rm -f ${gtflist}         

# Loop over samples
for sample in s1 s1t s2 s2t s3 s4 s5 s6 s7 s8 s9 s10 s11 s12

do
    # Define filepaths
    sample="${sample}.c2" # amend sample for all filepaths
    fastq="${data}/${sample}.fq"	# path to retrieve raw fastq 
    trim1="$data/${sample}.t1.fq"		# path to write/read adapter-trimmed fastq file 
    trim2="$data/${sample}.t2.fq"		# path to write/read  quality-trimmed fastq file 
    sam="${hisat_dir}/${sample}.sam"  # path to write/read  hisat2-generated SAM file
    bam="${hisat_dir}/${sample}.bam"  # path to write/read  samtools-converted BAM file 
    sorted_bam="${hisat_dir}/${sample}.sort.bam"	 # path to write/read  samtools-sorted BAM file    
    
    # Adapter trim with scythe
    scythe -o $trim1 -a illumina_adapter.fa -q sanger ${fastq}

    # Quality trim with sickle 
    sickle se -f $trim1 -t sanger -o $trim2  -q 10 -l 51 

    # Execute hisat2 to perform alignment
    hisat2 -p 4 --phred33 -x ${hs2index} -U ${trim2} -S ${sam}

    # Convert sam to bam with samtools
    samtools view -b -o ${bam} ${sam}
    # With samtools, sort bam
    samtools sort -o ${sorted_bam} ${bam}	 

    # Remove files that are no longer needed
    rm ${sam} ${bam}		 
    rm ${trim1} ${trim2}	

    # Create subdir to write stringtie results for this sample
    str_smp_dir="${stringtie_dir}/${sample}"		
    mkdir -p $str_smp_dir	

    sample_tr_gtf="${str_smp_dir}/${sample}_transcripts.gtf" # Filepath to output stringtie GTF
    # Assemble with stringtie
    stringtie -p 4 -t -e -B -G ${gtf} -o ${sample_tr_gtf} ${sorted_bam}	 

    # Append a line with the sample and stringtie GTF to gtflist
    gtfline="${sample} $sample_tr_gtf"
    echo ${gtfline} >> ${gtflist}

done 

# Run prepDE.py on the gtflist to get final trancript files
prepDE.py -i ${gtflist}
