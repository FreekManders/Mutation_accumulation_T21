#!/usr/bin/env python3

import sys
if sys.version_info[0] < 3:
    raise Exception("This script requires Python 3 to run.")

import argparse
import os

parser = argparse.ArgumentParser(description = "This script takes a somatic vcf file as its input (So filtered on non gt paramets). It tries to identify predisposition germline variants.")
parser.add_argument("-i", "--ini", required = True, help = "The ini file containing the settings of the script.")
args = parser.parse_args()

def ini_parser(ini_fname):
    r"""Parses a ini file with a key\tvalue format and returns a dictionary"""
    if not os.path.isfile(ini_fname): raise IOError("The ini file could not be read. This has caused an error.")
    ini_dict = {}
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
            ini_dict[key] = value
    return ini_dict

def istrue(flag):
    r"""Determine if a flag is set to true or not"""
    if flag.lower() in ["true", "t", "yes", "y", "do", "ok", "1"]:
        return True
    else:
        return False

def write_or_not(fname):
    write = not os.path.exists(fname) or overwrite
    return write

#Read in ini file
ini_dict = ini_parser(args.ini)
ini_needed = ["VCF", "OUT_DIR", "OVERWRITE", "sample", "GQ", "VAF", "DP", "snpsift", "R_predisposition_script", "R_predisposition_files_dir"]
missing = [i for i in ini_needed if i not in ini_dict.keys()]
if missing:
    raise ValueError("The following parameters are missing from the ini file: {0}".format(missing))


vcf = ini_dict["VCF"]
overwrite = istrue(ini_dict["OVERWRITE"])
dest = ini_dict["OUT_DIR"]
GQ = ini_dict["GQ"]
VAF = ini_dict["VAF"]
DP = ini_dict["DP"]
sample = ini_dict["sample"]
snpsift = ini_dict["snpsift"]
R_script = ini_dict["R_predisposition_script"]
R_predisp_files = ini_dict["R_predisposition_files_dir"]

#Finds location of sample. This is necessary, because snpsift crashes when told to filter on a sample whose name contains a "-"
with open(vcf) as vcf_file:
    for line in vcf_file:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            line = line.strip()
            header = line.split("\t")
            samples = header[9:]
            sample_loc = samples.index(sample)
            break

filtered_vcf = "{0}/{1}_filtered.vcf".format(dest, sample)
filter = "(GEN[{0}].GQ > {1}) & (GEN[{0}].DP > {2}) & ((GEN[{0}].AD[1]/(GEN[{0}].AD[0] + 0.001 + GEN[{0}].AD[1])) > {3})".format(sample_loc, GQ, DP, VAF)
if not os.path.isfile(filtered_vcf):
    print("Filtering on quality")
    os.system('java -Xmx4G -jar {0} filter "{1}" {2} > {3}'.format(snpsift, filter, vcf, filtered_vcf))

freq_vcf = "{0}/{1}_lowfreq.vcf".format(dest, sample)
if not os.path.isfile(freq_vcf):
    print("Filtering on frequency")
    os.system("""java -Xmx4G -jar {0} filter "(! dbNSFP_ExAC_AC > 10 ) & (! GoNLv5_AC > 10 )" {1} > {2}""".format(snpsift, filtered_vcf, freq_vcf))

impact_vcf = "{0}/{1}_highimpact.vcf".format(dest, sample)
if not os.path.isfile(impact_vcf):
    print("Filtering on impact")
    os.system("grep '#' {0} > {1}".format(freq_vcf, impact_vcf))
    os.system("grep -v '#' {0} | egrep 'HIGH|MODERATE' >> {1} ".format(freq_vcf, impact_vcf))

predisp_vcf = "{0}/{1}_predisposition.vcf".format(dest, sample)
if not os.path.isfile(predisp_vcf):
    print("Filtering on presence in gene lists")
    os.system("Rscript {0} --wdir {1} --sample {2} --pre_dir {3}".format(R_script, dest, sample, R_predisp_files))
