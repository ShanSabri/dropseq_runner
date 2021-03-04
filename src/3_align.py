#!/usr/bin/python

__author__  = "Shan Sabri"
__email__   = "ShanASabri@gmail.com"
__date__    = "09/08/2016"

import os

WRDIR = "/path/to/process/output/"
REFDIR = "/path/to/ref_data"

samples = [line.strip() for line in open(WRDIR + "/samples.txt")]

for idx, s in enumerate(samples):

    print "Processing: %s (%s/%s)" % (s, str(idx+1), str(len(samples)))

    TEMPDIR = WRDIR + "/" + s + "_temp/"
    OUTDIR  = WRDIR + "/" + s + "_output/"
    LOGDIR  = WRDIR + "/qsub_oe"
    INPUT = TEMPDIR + s + "_tagged_cell_umi_filtered_trimmed_smart_polyA_rmShort.fastq"
    bed = ""

    for GENOMEDIR in [REFDIR + "/Mus_musculus/UCSC/mm9", REFDIR + "/Homo_sapiens/UCSC/hg19"]:
        if not os.path.exists(TEMPDIR): os.makedirs(TEMPDIR)
        if not os.path.exists(OUTDIR): os.makedirs(OUTDIR)
        if not os.path.exists(LOGDIR): os.makedirs(LOGDIR)

        id = GENOMEDIR.split("/")[9]

	if id == "hg19": 
		bed = "/path/to/ref_data/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.bed"

	if id == "mm9":
		bed = "/path/to/ref_data/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/mm9.bed"

    ANNOTATION = WRDIR + "/" + id + "_longest_isoform.bed"
    OUTPUT_ALIGNED = OUTDIR + s + "_" + id + "_AlignedSorted"
    OUTPUT_BED = OUTDIR + s + "_" + id + "_AlignedSorted.bed"
    INTERSECT_OUTPUT = OUTDIR + s + "_" + id + "_Final.bed"
	DIST_OUT = OUTDIR + s + "_" + id + "_read_distribution.txt"

    align = "bowtie2 --time --threads 5 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50" + \
                     " -x " + GENOMEDIR + "/Sequence/Bowtie2Index/genome" + \
                     " -U " + INPUT + \
                     " -S " + TEMPDIR + s + "_" + id + "_" + "aligned.sam"

    bamtobed = "bedtools bamtobed -i " + OUTPUT_ALIGNED + ".bam > " + OUTPUT_BED

    map = "bedtools intersect -bed -u -wa -a " + OUTPUT_BED + " -b " + ANNOTATION + " > " + OUTDIR + id + "_temp.bed"
        
	join = "bedtools intersect -bed -wa -wb -a " + ANNOTATION + " -b " + OUTDIR + id + "_temp.bed > " + INTERSECT_OUTPUT
	
	intersect = "bedtools intersect -wa -wb -a " + OUTPUT_BED + " -b " + ANNOTATION + " > " + INTERSECT_OUTPUT

    os.system("echo '#!/bin/bash' > align_" + s + "_" + id)
    os.system("echo 'source ~/.bash_profile' >> align_" + s + "_" + id)
    os.system("echo '#$ -cwd' >> align_" + s + "_" + id)
    os.system("echo '#$ -o " + LOGDIR + "' >> align_" + s + "_" + id)
    os.system("echo '#$ -e " + LOGDIR + "' >> align_" + s + "_" + id)
    os.system("echo '#$ -m n' >> align_" + s + "_" + id)
    os.system("echo '#$ -l h_data=8G,h_rt=12:00:00' >> align_" + s + "_" + id)
    os.system("echo '#$ -pe shared 5' >> align_" + s + "_" + id)
    os.system("echo 'module load samtools' >> align_" + s + "_" + id)
	os.system("echo '" + align + "' >> align_" + s + "_" + id)
    os.system("echo 'samtools view -bS " + TEMPDIR + s + "_" + id + "_" + "aligned.sam | samtools sort - " + OUTPUT_ALIGNED + "' >> align_" + s + "_" + id)
	os.system("echo '" + bamtobed + "' >> align_" + s + "_" + id)
	os.system("echo '" + intersect + "' >> align_" + s + "_" + id) # correct
    os.system("qsub align_" + s + "_" +id)
    os.system("rm align_" + s + "_" +id)
