#!/usr/bin/python

__author__  = "Shan Sabri"
__email__   = "ShanASabri@gmail.com"
__date__    = "09/08/2016"

import datetime, editdistance, gzip 
from itertools import izip

wrdir        = "/path/to/qseqs/"
barcode_file = wrdir + "s_2_1.qseq.txt.gz"
idx_file     = wrdir + "s_2_2.qseq.txt.gz"
read_file    = wrdir + "s_2_3.qseq.txt.gz"

# Adjust idx as needed
true_idx     = [
    "TAAGGCGA", 
    "CGTACTAG", 
    "AGGCAGAA", 
    "TCCTGAGC", 
    "GGACTCCT", 
    "TAGGCATG"
]

out = {}
for t in true_idx:
    out[t] = (gzip.open(t + "_1.fastq.gz", "w"),
              gzip.open(t + "_2.fastq.gz", "w"))

start_time = datetime.datetime.now()

with gzip.open(barcode_file) as bf, gzip.open(idx_file) as idxf, gzip.open(read_file) as rf:

    print "START: ", start_time

    for lineNum, (b, i, r) in enumerate(izip(bf, idxf, rf)):

        b = b.strip().split("\t")
        i = i.strip().split("\t")
        r = r.strip().split("\t")

        obs_idx       = i[8]
        obs_read      = r[8]
        obs_read_qual = r[9]
        obs_bc        = b[8][0:20]
        obs_bc_qual   = b[9][0:20]
        read_seq_id   = "@" + ":".join(b[0:8]) + " length:" + str(len(obs_read))
        bc_seq_id     = "@" + ":".join(b[0:8]) + " length:" + str(len(obs_bc))

        def match_keys(_true_idx):
            return editdistance.eval(_true_idx, obs_idx)

        true_idx = min(out.keys(), key=match_keys)

        if editdistance.eval(true_idx, obs_idx) > 2:
            lineNum+=1
            continue

        print datetime.datetime.now()-start_time,"\tProcessing line", lineNum, ": ", ":".join(b[0:8]), obs_idx, " to ", true_idx

        out_barcode, out_read = out[true_idx]
        out_barcode.write("{0}\n{1}\n+\n{2}\n".format(bc_seq_id, obs_bc, obs_bc_qual))
        out_read.write("{0}\n{1}\n+\n{2}\n".format(read_seq_id, obs_read, obs_read_qual))

print "FINISH: ", datetime.datetime.now()
