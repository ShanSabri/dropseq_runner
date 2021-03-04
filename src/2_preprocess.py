#!/usr/bin/python

__author__  = "Shan Sabri"
__email__   = "ShanASabri@gmail.com"
__date__    = "09/08/2016"

import os

WRDIR = "/path/to/process/output/"
DST   = "/path/to/DST"

samples = [line.strip() for line in open(WRDIR + "/samples.txt")]

for idx, s in enumerate(samples):
    print "Processing: %s (%s/%s)" % (s, str(idx+1), str(len(samples)))

    TEMPDIR = WRDIR + "/" + s + "_temp/"
    OUTDIR  = WRDIR + "/" + s + "_output/"
    LOGDIR  = WRDIR + "/qsub_oe"

    if not os.path.exists(TEMPDIR): os.makedirs(TEMPDIR)
    if not os.path.exists(OUTDIR): os.makedirs(OUTDIR)
    if not os.path.exists(LOGDIR): os.makedirs(LOGDIR)

    tag_cells = DST + "/TagBamWithReadSequenceExtended " \
                      "INPUT=" + s + "_RAW_UNALIGNED.bam " \
                      "OUTPUT=" + TEMPDIR + s + "_tagged_cell.bam " \
                      "SUMMARY=" + OUTDIR + "/tagged_cell_summary.txt " \
                      "BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1"

    tag_molecules = DST + "/TagBamWithReadSequenceExtended " \
                          "INPUT=" + TEMPDIR + s + "_tagged_cell.bam " \
                          "OUTPUT=" + TEMPDIR + s + "_tagged_cell_umi.bam " \
                          "SUMMARY=" + OUTDIR + "/tagged_cell_umi_summary.txt " \
                          "BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1"

    filter_bam = DST + "/FilterBAM TAG_REJECT=XQ " \
                       "INPUT=" + TEMPDIR + s + "_tagged_cell_umi.bam " \
                       "OUTPUT=" + TEMPDIR + s + "_tagged_cell_umi_filtered.bam" \

    trim_starting_sequence = DST + "/TrimStartingSequence " \
                                   "INPUT=" + TEMPDIR + s + "_tagged_cell_umi_filtered.bam " \
                                   "OUTPUT=" + TEMPDIR + s + "_tagged_cell_umi_filtered_trimmed_smart.bam " \
                                   "OUTPUT_SUMMARY=" + OUTDIR + "/adapter_trimming_report.txt " \
                                   "SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5"

    trim_poly_a = DST + "/PolyATrimmer " \
                        "INPUT=" + TEMPDIR + s + "_tagged_cell_umi_filtered_trimmed_smart.bam " \
                        "OUTPUT=" + TEMPDIR + s + "_tagged_cell_umi_filtered_trimmed_smart_polyA.bam " \
                        "OUTPUT_SUMMARY=" + OUTDIR + "/polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=6"

    sam_to_fastq = "java -Xmx500m -jar /path/to/picard.jar SamToFastq " \
                        "INPUT=" + TEMPDIR + s + "_tagged_cell_umi_filtered_trimmed_smart_polyA.bam " \
                        "FASTQ=" + TEMPDIR + s + "_tagged_cell_umi_filtered_trimmed_smart_polyA.fastq " \

    os.system("echo '#!/bin/bash' > process_" + s)
    os.system("echo 'source ~/.bash_profile' >> process_" + s)
    os.system("echo '#$ -cwd' >> process_" + s)
    os.system("echo '#$ -o " + LOGDIR + "' >> process_" + s)
    os.system("echo '#$ -e " + LOGDIR + "' >> process_" + s)
    os.system("echo '#$ -m n' >> process_" + s)
    os.system("echo '#$ -l h_data=16G,h_rt=12:00:00' >> process_" + s)
    os.system("echo '" + tag_cells + "' >> process_" + s)
    os.system("echo '" + tag_molecules + "' >> process_" + s)
    os.system("echo '" + filter_bam + "' >> process_" + s)
    os.system("echo '" + trim_starting_sequence + "' >> process_" + s)
    os.system("echo '" + trim_poly_a + "' >> process_" + s)
    os.system("echo '" + sam_to_fastq + "' >> process_" + s)
    os.system("qsub process_" + s)
    os.system("rm process_" + s)
