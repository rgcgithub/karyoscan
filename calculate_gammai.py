#!/usr/bin/python

import os
import sys
import math
import concurrent.futures
import multiprocessing
from itertools import izip

usage_str = "<windows_file> <mappability_threshold [0-1]> <min_gc [0-1]> <max_gc [0-1]> <sample_name:coverage_file pairs input TSV file> <output_filename>"

if len(sys.argv) < 5:
    print("Usage: %s %s" % (sys.argv[0], usage_str))
    sys.exit(0)

windows_path = sys.argv[1]
mappability = float(sys.argv[2])
min_gc = float(sys.argv[3])
max_gc = float(sys.argv[4])
sample_list = sys.argv[5]
outfile = sys.argv[6]

sample_names = []
sample_infiles = []

with open(sample_list, 'r') as fin:
    for line in fin:
        sid, covfile = [line.rstrip().split('\t')[i] for i in (0,1)]
        sample_names.append(sid)
        sample_infiles.append(covfile)
fin.close()

def process_one(sample_id, coverage_file, windows):
    ri_values = {"1":0,"2":0,"3":0,"4":0,"5":0,"6":0,"7":0,"X":0,"8":0,"9":0,"10":0,"11":0,"12":0,"13":0,"14":0,"15":0,"16":0,"17":0,"18":0,"20":0,"Y":0,"19":0,"22":0,"21":0}

    with open(windows, 'r') as fwindows, \
         open(coverage_file, 'r') as fcov:
         
        for i, j in izip(fwindows, fcov):
            window_tokens = i.split('\t')
            coverage_tokens = j.split('\t')
            if not window_tokens[0:3] == coverage_tokens[0:3]:
                print("Sample %s failed due to mismatching BED records between the coverage file (%s: %s) and windows file (%s: %s)!" % (sample_id, coverage_file, ":".join(coverage_tokens[0:3]), windows, ":".join(window_tokens[0:3])))
                return(-1)
            coverage = float(coverage_tokens[3])
            if not math.isnan(coverage) \
                and int(window_tokens[-1]) == 3 \
                and float(window_tokens[-2]) >= mappability \
                and float(window_tokens[-3]) >= min_gc \
                and float(window_tokens[-3]) <= max_gc:
                    ri_values[str(window_tokens[0])] += coverage
    fcov.close()
    fwindows.close()
 
    autosome_sum = float(sum(ri_values.values()) - ri_values['X'] - ri_values['Y'])
    gammai_values = {k: float(v) / float(autosome_sum) if k == "X" or k == "Y" else float(v) / float(autosome_sum-v) for k,v in ri_values.items()}
   
    return (sample_id, "\t".join("{:.8g}".format(val) for key,val in sorted(gammai_values.items())))


nproc = multiprocessing.cpu_count()

output_str = ""
with concurrent.futures.ThreadPoolExecutor(max_workers=4*nproc) as proc_executor:
    proc_futs = [proc_executor.submit(process_one, sample_names[j], sample_infiles[j], windows_path) for j in xrange(len(sample_names))]
    for pfut in concurrent.futures.as_completed(proc_futs):
        if not pfut.result() == -1: # -1 = Sample failed
            output = pfut.result()
            sample_id = output[0]
            output_str = output_str + "%s\t%s\n" % (sample_id, output[1])

concurrent.futures.wait(proc_futs, return_when=concurrent.futures.ALL_COMPLETED)

fout = open(outfile, 'w')
fout.write("SampleID\tchr" + "\tchr".join(sorted(str(i) for i in range(1,23))) + "\tchrX\tchrY\n")
fout.write(output_str)
