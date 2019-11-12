#Creates a shell script submitting telomerecat for each bam file

import os

#Data
bam_dirs = ["/hpc/pmc_vanboxtel/processed/fetal_clones/OS1/bams/",
"/hpc/pmc_vanboxtel/processed/fetal_clones/NR1/bams/",
"/hpc/pmc_vanboxtel/processed/fetal_clones/NR2/bams/",
"/hpc/pmc_vanboxtel/processed/fetal_clones/N01/bams/",
"/hpc/pmc_vanboxtel/processed/fetal_clones/MH3/bams/",
"/hpc/pmc_vanboxtel/processed/HMFreg0367/MH2/MH2LIMPPCD13/mapping/",
"/hpc/pmc_vanboxtel/processed/HMFreg0367/MH2/MH2LIMPPCL13/mapping/",
"/hpc/pmc_vanboxtel/processed/HMFreg0367/MH2/MH2SIBULK/mapping/"]

wd = "/hpc/pmc_vanboxtel/projects/Freek_trees/telomeres/telomerecat/"

#Set wd
os.chdir(wd)

#Determine file paths
files = []
for dir in bam_dirs:
    for r, d, f in os.walk(dir):
        for file in f:
            if file.endswith(".bam"):
                files.append(os.path.join(r, file))

#Create shell script and run
for bam_file in files:
    base = os.path.basename(bam_file)
    base = base.replace(".bam", "")
    sample = base.replace("_dedup.realigned", "")
    bash_text = """#!/bin/bash\n\n
#$ -N Telomerecat\n
#$ -cwd\n#$ -pe threaded 8\n
#$ -l h_rt=2:0:0\n
#$ -o {0}_telomercat_output.txt\n
#$ -e {0}_telomercat_error.txt\n\n
cd {1}\n
. /hpc/pmc_vanboxtel/tools/telomerecat/venv_2.7.12/bin/activate\n\n
BAM={2}\n\n
echo $BAM\n
telomerecat bam2length -p8 $BAM\n
echo Done""".format(sample, wd, bam_file)
    bash_fname = "qsub_telomerecat_" + sample + ".sh"
    if not os.path.isfile(bash_fname):#Do not run, when it already exists.
        with open(bash_fname, "w") as bash:
            bash.write(bash_text)
        os.system("qsub {0}".format(bash_fname))
