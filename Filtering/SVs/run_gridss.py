#!/usr/bin/env python3

import sys
if sys.version_info[0] < 3:
    raise Exception("This script requires Python 3 to run.")

import os
from os.path import join as j
import argparse
from collections import defaultdict
from shutil import copy

parser = argparse.ArgumentParser(description = "This script filters a somatic vcf created by the hmf-pipeline. It can also be used to filter a vcf from the iap. It only filters on number of alleles, blacklists and mq.")
parser.add_argument("-i", "--ini", required = True, help = "The ini file containing the settings of the script.")
args = parser.parse_args()

def ini_parser(ini_fname):
    r"""Parses a ini file with a key\tvalue format and returns a dictionary"""
    if not os.path.isfile(ini_fname): raise IOError("The ini file does not exist. This has caused an error.")
    ini_dict = defaultdict(list)
    with open(ini_fname) as ini:
        for line in ini:
            line = line.strip()
            if line.startswith("#") or line == "":
                continue
            if "\t" in line:
                key, value = line.split("\t")
            else:
                key = line
                value = ""
            if key == "BAM":
                ini_dict[key].append(value)
            else:
                ini_dict[key] = value
    return ini_dict



#Read in ini file
ini_dict = ini_parser(args.ini)
ini_needed = ["BLACKLIST", "REFERENCE", "OUTDIR", "TIME", "RAM", "THREADS", "EMAIL", "SAMPLE_NAME", "Guix_profile", "GRIDSS_JAR", "VAF", "REF_CELL", "PON_DIR", "Rscript", "cancer_genes"]
missing = [i for i in ini_needed if i not in ini_dict.keys()]
if missing:
    raise ValueError("The following parameters are missing from the ini file: {0}".format(missing))

#Create output directories
outdir = ini_dict["OUTDIR"]
if not os.path.isdir(outdir):
    os.mkdir(outdir)

#Copy the ini to the output folder and set the output as the working directory.
copy(args.ini, outdir)
os.chdir(outdir)
#tmp_dir = j(outdir, "tmp")
#if not os.path.isdir(tmp_dir):
#    os.mkdir(tmp_dir)


#Get guix profile and necessary jars.
Guix_profile = ini_dict["Guix_profile"]
Guix_profile_source = j(Guix_profile, "etc", "profile") #Used to source the guix profile
#gridss = j(Guix_profile, "share", "java", "gridss", "gridss.jar")
gridss = ini_dict["GRIDSS_JAR"]
snpsift = j(Guix_profile, "share", "java", "snpeff", "SnpSift.jar")

#Determine name of main output files
sample = ini_dict["SAMPLE_NAME"]
output_sample = j(outdir, sample)
output = output_sample + ".sv.vcf"
assembly = output_sample + ".gridss.assembly.bam"

#Create the gridss command
ram = ini_dict["RAM"]
threads = ini_dict["THREADS"]
total_ram = int(ram) * int(threads)
ram_gridss = int(ram) - 10
reference = ini_dict["REFERENCE"]
blacklist = ini_dict["BLACKLIST"]
email = ini_dict["EMAIL"]
time = ini_dict["TIME"]

part1 = """java -ea -Xmx{0}G \
-Dsamjdk.create_index=true \
-Dsamjdk.use_async_io_read_samtools=true \
-Djava.io.tmpdir=$TMPDIR \
-Dsamjdk.use_async_io_write_samtools=true \
-Dsamjdk.use_async_io_write_tribble=true \
-Dsamjdk.buffer_size=4194304 \
-Dgridss.gridss.output_to_temp_file=true \
-cp {1} gridss.CallVariants \
TMP_DIR=$TMPDIR \
WORKING_DIR={2} \
REFERENCE_SEQUENCE="{3}" \
OUTPUT="{4}" \
ASSEMBLY="{5}" \
BLACKLIST="{6}" \
WORKER_THREADS={7}""".format(ram_gridss, gridss, outdir, reference, output, assembly, blacklist, threads)


if "BAMS" in ini_dict.keys():
    bam_dir = ini_dict["BAMS"]
    f_list = []
    for (dirpath, dirnames, filenames) in os.walk(bam_dir):
        f_list.extend(filenames)
        break
    f_list = [f for f in f_list if f.endswith(".bam")]
    bams = [j(bam_dir, f) for f in f_list]
else:
    bams = ini_dict["BAM"]

for bam in bams:
    part1 = part1 + " INPUT=" + bam

    bam_basename = os.path.basename(bam)
    bam_basename = bam_basename[:-4]
    bam_name = bam_basename.split("_", 1)[0]
    part1 = part1 + " INPUT_LABEL=" + bam_name

gridss_command = part1

#Create the snpsift filter command
filtered_vcf = output_sample + ".sv.filtered.vcf"
filter = "(FILTER = 'PASS') & (QUAL >= 1000) & (AS > 0) & (RAS > 0)"
snpsift_command = 'java -Xmx{0}G -jar {1} filter "{2}" {3} > {4}'.format(ram_gridss, snpsift, filter, output, filtered_vcf)

#Create the R filtering command
R_command = "Rscript {0} --SAMPLE_NAME {1} --DIR {2} --pon_dir {3} --vaf {4} --ref_cel {5} --cancer_genes {6}".format(ini_dict["Rscript"], sample, outdir, ini_dict["PON_DIR"], ini_dict["VAF"], ini_dict["REF_CELL"], ini_dict["cancer_genes"])

#Create a qsub submission script
gridss_bash_fname = j(outdir, "run_gridss.sh")
out_log_fname = j(outdir, "run_gridss_output.txt")
error_log_fname = j(outdir, "run_gridss_error.txt")
jobname = sample + "_rungridss"
gridss_bash = "#!/bin/bash\n\n#$ -S /bin/bash\n#$ -l h_rt={0}\n#$ -l h_vmem={1}G\n#$ -pe threaded {2}\n#$ -l tmpspace=200G\n#$ -cwd\n#$ -o {3}\n#$ -e {4}\n#$ -M {5}\n#$ -m eas\n#$ -N {6}\n\nsource {7}\n\n{8}\n\n{9}\n\n{10}".format(time, ram, threads, out_log_fname, error_log_fname, email, jobname, Guix_profile_source, gridss_command, snpsift_command, R_command)
with open(gridss_bash_fname, "w") as gridss_bash_file:
    gridss_bash_file.write(gridss_bash)
os.system("qsub {0}".format(gridss_bash_fname))
print("Submitted gridss job for sample: {0}".format(sample))
