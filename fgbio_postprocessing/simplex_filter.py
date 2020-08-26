#!/home/arorak/miniconda2/bin/python2.7

import re
import sys
import pysam
import os


USAGE = "USAGE: " + sys.argv[0] + " <INPUT_BAM_FILE> [ <OUTPUT_BAM_FILE> ]\n"
min_depth = 3

if len(sys.argv) < 2:
    sys.stderr.write(USAGE)
    sys.exit(1)

inputbam = sys.argv[1]

if not os.path.isfile(inputbam):
    sys.stderr.write("Input BAM file {0} does not exist.\n".format(inputbam))
    sys.exit(1)

if len(sys.argv) > 2:
    simplexbam = sys.argv[2]
else:
    simplexbam =  re.sub(".bam$", "", inputbam) + ".simplex.bam" #making sure suffix is added even if the input file doesn't end with .bam suffix.

bamfile = pysam.AlignmentFile(inputbam, "rb")
simplex = pysam.AlignmentFile(simplexbam, "wb", template=bamfile)

for read in bamfile.fetch():
    try:
        strand1_dp=read.get_tag("aD")
        strand2_dp=read.get_tag("bD")
        consens_dp=read.get_tag("cD")
        # check duplex conditions
        if consens_dp >= min_depth and min(strand1_dp, strand2_dp) == 0:
            simplex.write(read)
    except:
        continue

bamfile.close()
simplex.close()
try:
    pysam.index(simplexbam, re.sub(".bam$", "", simplexbam) + ".bai")
except:
    sys.stderr.write("Could not index Simplex bam file {0}.\n".format(simplexbam))
