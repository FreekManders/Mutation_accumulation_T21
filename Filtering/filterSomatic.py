#!/usr/bin/env python3

import sys
if sys.version_info[0] < 3:
    raise Exception("This script requires Python 3 to run.")

import argparse
import os

parser = argparse.ArgumentParser(description = "This script filters a somatic vcf created by the hmf-pipeline. It can also be used to filter a vcf from the iap. It only filters on number of alleles, blacklists and mq.")
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
ini_needed = ["FILE", "OUT_DIR", "OVERWRITE", "only_snv", "only_indel", "chroms", "sample_name", "qual"]
missing = [i for i in ini_needed if i not in ini_dict.keys()]
if missing:
    raise ValueError("The following parameters are missing from the ini file: {0}".format(missing))

#Filter on quality
vcf = ini_dict["FILE"]
overwrite = istrue(ini_dict["OVERWRITE"])

if not os.path.isdir(ini_dict["OUT_DIR"]):
    os.mkdir(ini_dict["OUT_DIR"])

if vcf.endswith(".gz"):
    input = "--gzvcf {0}".format(vcf)
else:
    input = "--vcf {0}".format(vcf)

extra_filters = ""

if ini_dict["chroms"] is not "":
    chroms = ini_dict["chroms"].split(",")
    for chrom in chroms:
        extra_filters = extra_filters + "--chr {0} ".format(chrom)

if istrue(ini_dict["only_snv"]):
    extra_filters = extra_filters + "--remove-indels "

if istrue(ini_dict["only_indel"]):
    extra_filters = extra_filters + "--keep-only-indels "

if ini_dict["qual"].lower() not in ["false", "f", "0", ""]:
    extra_filters = extra_filters + "--minQ " + ini_dict["qual"] + " "


if ini_dict["max_alleles"].lower() not in ["false", "f", "0", ""]:
    extra_filters = extra_filters + "--max-alleles " + ini_dict["max_alleles"] + " "

out_name = ini_dict["sample_name"] + "_somatic_filtered"
out_vcf = os.path.join(ini_dict["OUT_DIR"], out_name)
out_vcf_fullname = out_vcf + ".recode.vcf"
if write_or_not(out_vcf_fullname):
    command = "vcftools {0} --remove-filtered-all {1}--recode --recode-INFO-all --out {2}".format(input, extra_filters, out_vcf)
    print(command)
    os.system(command)
    print("Filtered the somatic vcf for sample: {0} on quality".format(ini_dict["sample_name"]))
else:
    print("Quality filtering on the somatic vcf for sample: {0} already performed".format(ini_dict["sample_name"]))


#Filter on blacklist
in_vcf = out_vcf_fullname
out_name = ini_dict["sample_name"] + "_somatic_filtered_noblacklist.vcf"
nonblacklist_vcf = os.path.join(ini_dict["OUT_DIR"], out_name)
if ini_dict["blacklists"] is not "" and write_or_not(nonblacklist_vcf):
    blacklists = ini_dict["blacklists"].split(",")
    i = 1
    for blacklist in blacklists:
        tmp_out = os.path.join(ini_dict["OUT_DIR"], "temp{0}".format(i))
        tmp_out_fullname = tmp_out + ".recode.vcf"
        command = "vcftools --vcf {0} --exclude-positions {1} --recode --recode-INFO-all --out {2}".format(in_vcf, blacklist, tmp_out)
        os.system(command)
        in_vcf = tmp_out_fullname
        i = i + 1
    #Change name of final vcf
    final_vcf_oldname = tmp_out + ".recode.vcf" #This should be the same as the final in_vcf, but in case of errors this won't override the original in_vcf
    os.system("mv {0} {1}".format(final_vcf_oldname, nonblacklist_vcf))
    print("Filtered the somatic vcf on the blacklists for sample: {0}".format(ini_dict["sample_name"]))
    files = os.listdir(ini_dict["OUT_DIR"])
    files_path = [os.path.join(ini_dict["OUT_DIR"], file) for file in files]
    nothing = [os.remove(i) for i in files_path if "temp" in i]
else:
    print("Either no blacklist was given or sample: {0} was already filtered on blacklists".format(ini_dict["sample_name"]))

#Filter on mq
out_name = ini_dict["sample_name"] + "_somatic_filtered_noblacklist_MQ.vcf"
mqfiltered_vcf = os.path.join(ini_dict["OUT_DIR"], out_name)
if ini_dict["MQ"].lower() not in ["false", "f", "0", ""] and write_or_not(mqfiltered_vcf):
    filter = "(MQ >= {0})".format(ini_dict["MQ"])
    command = 'java -Xmx4G -jar {0} filter "{1}" {2} > {3}'.format(ini_dict["snpsift"], filter, nonblacklist_vcf, mqfiltered_vcf)
    os.system(command)
    print("Filtered the somatic vcf on MQ for sample: {0}".format(ini_dict["sample_name"]))
else:
    print("MQ filtering on the somatic vcf for sample: {0} already performed".format(ini_dict["sample_name"]))
