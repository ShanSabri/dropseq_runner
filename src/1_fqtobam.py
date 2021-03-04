#!/usr/bin/python

__author__  = "Shan Sabri"
__email__   = "ShanASabri@gmail.com"
__date__    = "09/08/2016"

import os

FQDIR      = "/path/to/demux/output/"
WRDIR      = "/path/to/process/output/"
LOGDIR     = WRDIR + "/qsub_oe/"
FASTQTOSAM = "java -jar /path/to/picard.jar FastqToSam"

if not os.path.exists(LOGDIR): os.makedirs(LOGDIR)

samples = [line.strip() for line in open(WRDIR + "/samples.txt")]

for idx, s in enumerate(samples):
    os.system("printf 'Processing: %s (%s/%s)\n' " + s + " " + str(idx+1) + " " + str(len(samples)))
    os.system("echo '#!/bin/csh' > runPicardFastqToSam_" + s)
    os.system("echo 'source ~/.bash_profile' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -cwd' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -o " + LOGDIR + "' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -e " + LOGDIR + "' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -m n' >> runPicardFastqToSam_" + s)
    os.system("echo '#$ -l h_data=16G,h_rt=8:00:00' >> runPicardFastqToSam_" + s)
    os.system("echo '" + FASTQTOSAM + " FASTQ="  + FQDIR + "/" + s + "_1.fq.gz" +
                                      " FASTQ2=" + FQDIR + "/" + s + "_2.fq.gz" +
                                      " OUTPUT=" + WRDIR + "/" + s + "_RAW_UNALIGNED.bam"     +
                                      " SAMPLE_NAME=" + s + "' >> runPicardFastqToSam_" + s)
    os.system("qsub runPicardFastqToSam_" + s)
    os.system("rm runPicardFastqToSam_" + s)
