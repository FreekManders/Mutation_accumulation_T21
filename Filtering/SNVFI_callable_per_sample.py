#Run SNVFI filtering, germline filtering, CallableLoci filtering and more for all chosen samples. Filtersettings are changed based on the sex and trisomy status of samples.

#Install dependencies with: guixr package -i r@3.5.1 r-genomicranges@1.34.0 grep@3.1 coreutils@8.30 r-tidyverse@1.2.1 r-variantannotation@1.28.10 r-biocinstaller@1.32.1 r-ggplot2@3.1.0 r-reshape2@1.4.3 r-optparse@1.6.1 vcftools@0.1.15 snpeff@4.3t bedtools@2.27.1 bio-vcf@0.9.2 python@3.7.0 r-gdata@2.18.0
#To run: source ~/guix_SNVFI_callable_per_sample/etc/profile
# python3.7 /hpc/pmc_vanboxtel/tools/scripts/Freek_misc/SNVFI_callable_per_sample.py

import sys
if sys.version_info[0] < 3:
    raise Exception("This script requires Python 3 to run.")

import os

out = "/hpc/pmc_vanboxtel/projects/Freek_subclones_fetalt21/"
snvfi = "/hpc/pmc_vanboxtel/tools/SNVFI/SNVFI_run_v1.2.sh"
snvfi_config = "/hpc/pmc_vanboxtel/tools/SNVFI/SNVFI_default_Freek_noSGE.config"
indelfi = "/hpc/pmc_vanboxtel/tools/INDELFI/INDELFI_v1.4.0.pl"
snpsift = "/home/pmc_research/fmanders/Freek_guix-profile/share/java/snpeff/SnpSift.jar"
R_predisposition_script = "/hpc/pmc_vanboxtel/tools/scripts/Freek_misc/predisposistion_genes2.R"
R_predisposition_files_dir = "/hpc/pmc_vanboxtel/data/Predisposition/"
guix_profile_source = "/home/pmc_research/fmanders/guix_SNVFI_callable_per_sample/etc/profile"

#fetuses_dir = [{"dir_in" : "/hpc/pmc_vanboxtel/processed/Fetal_Body_Map/F100916W15/", "fetus" : "F100916W15", "bulk" : "F100916W15-REF", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/Fetal_Body_Map/F100916W17/", "fetus" : "F100916W17", "bulk" : "F100916W17-REF", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/Fetal_Body_Map/E080416/", "fetus" : "E080416", "bulk" : "E080416SKIN", "sex" : "M"}]
#fetuses_dir = [{"dir_in" : "/hpc/pmc_vanboxtel/processed/HMFreg0349/N01/", "fetus" : "N01", "bulk" : "N01SKIN", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/HMFreg0349/NR1/", "fetus" : "NR1", "bulk" : "NR1SKIN", "sex" : "M"}]
#fetuses_dir = [{"dir_in" : "/hpc/pmc_vanboxtel/processed/HMFreg0349/N01/", "fetus" : "N01", "bulk" : "N01SKIN", "sex" : "M"}]
#fetuses_dir = [{"dir_in" : "/hpc/pmc_vanboxtel/processed/HMFreg0349/N01/", "fetus" : "N01", "bulk" : "N01SKIN", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/HMFreg0349/NR1/", "fetus" : "NR1", "bulk" : "NR1SKIN", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/HMFreg0335/NC1/", "fetus" : "NC1", "bulk" : "NC1SIBULK", "sex" : "F"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/HMFreg0335/NR2/", "fetus" : "NR2", "bulk" : "NR2SKIN", "sex" : "M"}]
#fetuses_dir = [{"dir_in" : "/hpc/pmc_vanboxtel/processed/fetal_clones/N01/output/", "fetus" : "N01", "bulk" : "N01skin_combined", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/fetal_clones/NR1/output/", "fetus" : "NR1", "bulk" : "NR1skin_combined", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/fetal_clones/NR2/output/", "fetus" : "NR2", "bulk" : "NR2skin_combined", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/HMFreg0367/MH2/", "fetus" : "MH2", "bulk" : "MH2SIBULK", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/HMFreg0367/MH3/", "fetus" : "MH3", "bulk" : "MH3VILI", "sex" : "F"}]
#fetuses_dir = [{"dir_in" : "/hpc/pmc_vanboxtel/processed/HMFreg0409/OS1/", "fetus" : "OS1", "bulk" : "OS1SKINREF", "sex" : "F"}]
#fetuses_dir = [{"dir_in" : "/hpc/pmc_vanboxtel/processed/HMFreg0409/TMD14438/", "fetus" : "TMD14438", "bulk" : "TMD14438TCELLBULK", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/fetal_tmd/TMD14743/output/", "fetus" : "TMD14743", "bulk" : "TMD14743BCELL", "sex" : "F"}]
#fetuses_dir = [{"dir_in" : "/hpc/pmc_vanboxtel/processed/fetal_clones/MH3/output/", "fetus" : "MH3", "bulk" : "MH3SIBULKREF", "sex" : "F"}]
#fetuses_dir = [{"dir_in" : "/hpc/pmc_vanboxtel/processed/fetal_clones/N01/output/", "fetus" : "N01", "bulk" : "N01skin_combined", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/fetal_clones/NR1/output/", "fetus" : "NR1", "bulk" : "NR1skin_combined", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/fetal_clones/NR2/output/", "fetus" : "NR2", "bulk" : "NR2skin_combined", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/HMFreg0367/MH2/", "fetus" : "MH2", "bulk" : "MH2SIBULK", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/fetal_clones/MH3/output/", "fetus" : "MH3", "bulk" : "MH3SIBULKREF", "sex" : "F"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/fetal_clones/OS1/output/", "fetus" : "OS1", "bulk" : "OS1SKINREF", "sex" : "F"}]
fetuses_dir = [{"dir_in" : "/hpc/pmc_vanboxtel/processed/fetal_subclones/N01/output/", "fetus" : "N01", "bulk" : "N01skin_combined", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/fetal_subclones/NR1/output/", "fetus" : "NR1", "bulk" : "NR1skin_combined", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/fetal_subclones/NR2/output/", "fetus" : "NR2", "bulk" : "NR2skin_combined", "sex" : "M"}, {"dir_in" : "/hpc/pmc_vanboxtel/processed/fetal_subclones/OS1/output/", "fetus" : "OS1", "bulk" : "OS1SKINREF", "sex" : "F"}]




#loop over fetuses
for dir in fetuses_dir:
    dir_in = dir["dir_in"]
    fetus = dir["fetus"]
    bulk = dir["bulk"]
    sex = dir["sex"]

    #Create output directory for fetus
    out_dir = out + fetus + "/"
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    #Identify the vcf
    vcf = dir_in + fetus + ".filtered_variants_snpEff_snpSift_Cosmicv76_GoNLv5.vcf"
    vcf_copied = out + fetus + "/" + fetus + ".filtered_variants_snpEff_snpSift_Cosmicv76_GoNLv5.vcf"
    if not os.path.isfile(vcf):
        vcf = dir_in + "output.filtered_variants_snpEff_snpSift_Cosmicv76_GoNLv5.vcf"
        vcf_copied = out + fetus + "/output.filtered_variants_snpEff_snpSift_Cosmicv76_GoNLv5.vcf"
    if not os.path.isfile(vcf) and os.path.isfile(vcf_copied): #If vcf in processed is removed, use the copied file
        vcf = vcf_copied

    #Copy the vcf to the output directory. This way the vcf can still be used if the original iap output has been removed.
    if not os.path.isfile(vcf_copied):
        print("Copying the vcf for fetus {0} to the output directory.".format(fetus))
        os.system("cp {0} {1}".format(vcf, vcf_copied))
        print("Copied the vcf for fetus {0} to the output directory.".format(fetus))


    #Create a list of the sample names in the vcf. This is used to determine which column contains each sample / bulk
    with open(vcf) as vcf_file:
        for line in vcf_file:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                line = line.strip()
                header = line.split("\t")
                samples_bulk = header[9:]
                samples = samples_bulk.copy()
                samples.remove(bulk)
                break


    #Create sub output directory for callableloci
    out_callable = out_dir + "callableloci/"
    if not os.path.isdir(out_callable):
        os.mkdir(out_callable)

    #Get callable loci for the bulk
    bulk_bed = dir_in + bulk + "/" + bulk + "_CallableLoci.bed"
    bulk_callable_fname = out_callable + bulk + "CallableLoci_CALLABLE.bed"
    if not os.path.isfile(bulk_bed):
        bulk_bed = dir_in + bulk + "_dedup.realigned/" + bulk + "_dedup.realigned_CallableLoci.bed"
    if os.path.isfile(bulk_bed) and not os.path.isfile(bulk_callable_fname):
        command = "grep 'CALLABLE' {0} > {1}".format(bulk_bed, bulk_callable_fname)
        os.system(command)
        print("Filtered CALLABLE bulk file for fetus: {0}".format(fetus))

    bulk_autosomal_fname = out_callable + bulk + "CallableLoci_autosomal.bed"
    if os.path.isfile(bulk_callable_fname) and not os.path.isfile(bulk_autosomal_fname):
        command = "sed '/[XYMT]/d' {0} > {1}".format(bulk_callable_fname, bulk_autosomal_fname)
        os.system(command)
        print("Filtered for autosomal regions on the CALLABLE bulk file for fetus: {0}".format(fetus))

    bulk_autosomalx_fname = out_callable + bulk + "CallableLoci_autosomal_X.bed"
    if os.path.isfile(bulk_callable_fname) and not os.path.isfile(bulk_autosomalx_fname):
        command = "sed '/[YMT]/d' {0} > {1}".format(bulk_callable_fname, bulk_autosomalx_fname)
        os.system(command)
        print("Filtered for autosomal and x regions on the CALLABLE bulk file for fetus: {0}".format(fetus))

    #Run filterSomatic.py, whose output will later be used for finding shared mutations
    out_dir_filtersomatic = out_dir + "shared"
    if not os.path.isdir(out_dir_filtersomatic):
        os.mkdir(out_dir_filtersomatic)
    ini_fname = out_dir + "filterSomatic.ini"
    if not os.path.isfile(ini_fname):
        ini_text = "### .ini file for the filterSomatic.py description\n\nFILE\t{0}\nOUT_DIR\t{1}\nOVERWRITE\tT\nonly_snv\tT\nonly_indel\tF\nqual\tF\nmax_alleles\t2\nchroms\t1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X\nsample_name\t{2}\nblacklists\t/hpc/pmc_vanboxtel/data/Mutation_blacklists/MSC_healthyBM_raw_variants.vcf,/hpc/pmc_vanboxtel/data/dbSNP/nocosmic_dbsnp_137.b37.vcf\n\nsnpsift\t/home/pmc_research/fmanders/Freek_guix-profile/share/java/snpeff/SnpSift.jar\nMQ\t60".format(vcf, out_dir_filtersomatic, fetus)
        with open(ini_fname, "w") as ini:
            ini.write(ini_text)
        command = "python3.7 /hpc/pmc_vanboxtel/tools/filterSomatic/filterSomatic.py -i {0}".format(ini_fname)
        os.system(command)

    # #Look for germline predisposition mutations in the filtered vcf.
    # somatic_vcf_fname = out_dir_filtersomatic + "/" + fetus + "_somatic_filtered_noblacklist_MQ.vcf"
    # out_dir_germline = out_dir + "germline"
    # if not os.path.isdir(out_dir_germline):
    #     os.mkdir(out_dir_germline)
    # ini_fname = out_dir_germline + "/Filtergermline_impact.ini"
    # if not os.path.isfile(ini_fname):
    #     ini_text = "### .ini file for Filtergermline_impact2.py\n\nVCF\t{0}\nOUT_DIR\t{1}\nOVERWRITE\tT\nsample\t{2}\nGQ\t50\nVAF\t0.3\nDP\t10\n\nsnpsift\t{3}\nR_predisposition_script\t{4}\nR_predisposition_files_dir\t{5}".format(somatic_vcf_fname, out_dir_germline, bulk, snpsift, R_predisposition_script, R_predisposition_files_dir)
    #     with open(ini_fname, "w") as ini:
    #         ini.write(ini_text)
    #     command = "python3.7 /hpc/pmc_vanboxtel/tools/scripts/Freek_misc/Filtergermline_impact2.py -i {0}".format(ini_fname)
    #     os.system(command)

    #Run filterSomatic.py for indels, whose output will later be used for finding shared mutations
    out_dir_filtersomatic = out_dir + "indels"
    if not os.path.isdir(out_dir_filtersomatic):
        os.mkdir(out_dir_filtersomatic)
    ini_fname = out_dir + "filterSomatic_indels.ini"
    if not os.path.isfile(ini_fname):
        ini_text = "### .ini file for the filterSomatic.py description\n\nFILE\t{0}\nOUT_DIR\t{1}\nOVERWRITE\tT\nonly_snv\tF\nonly_indel\tT\nqual\t250\nmax_alleles\t2\nchroms\t1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X\nsample_name\t{2}\nblacklists\t/hpc/pmc_vanboxtel/data/Mutation_blacklists/MSC_healthyBM_raw_variants.vcf,/hpc/pmc_vanboxtel/data/dbSNP/nocosmic_dbsnp_137.b37.vcf\n\nsnpsift\t/home/pmc_research/fmanders/Freek_guix-profile/share/java/snpeff/SnpSift.jar\nMQ\t60".format(vcf, out_dir_filtersomatic, fetus)
        with open(ini_fname, "w") as ini:
            ini.write(ini_text)
        command = "python3.7 /hpc/pmc_vanboxtel/tools/filterSomatic/filterSomatic.py -i {0}".format(ini_fname)
        os.system(command)

    # #Look for germline predisposition indels in the filtered vcf.
    # somatic_vcf_fname = out_dir_filtersomatic + "/" + fetus + "_somatic_filtered_noblacklist_MQ.vcf"
    # out_dir_germline = out_dir + "germline_indel"
    # if not os.path.isdir(out_dir_germline):
    #     os.mkdir(out_dir_germline)
    # ini_fname = out_dir_germline + "/Filtergermline_impact.ini"
    # if not os.path.isfile(ini_fname):
    #     ini_text = "### .ini file for Filtergermline_impact2.py\n\nVCF\t{0}\nOUT_DIR\t{1}\nOVERWRITE\tT\nsample\t{2}\nGQ\t50\nVAF\t0.3\nDP\t10\n\nsnpsift\t{3}\nR_predisposition_script\t{4}\nR_predisposition_files_dir\t{5}".format(somatic_vcf_fname, out_dir_germline, bulk, snpsift, R_predisposition_script, R_predisposition_files_dir)
    #     with open(ini_fname, "w") as ini:
    #         ini.write(ini_text)
    #     command = "python3.7 /hpc/pmc_vanboxtel/tools/scripts/Freek_misc/Filtergermline_impact2.py -i {0}".format(ini_fname)
    #     os.system(command)


#Run SNVFI per sample
    samples = [sample for sample in samples if not sample[-1].isdigit()] #Run only for subclones with a letter at the end of the name
    for sample in samples:
        out_dir_snvfi = out_dir + sample + "/"
        if not os.path.isdir(out_dir_snvfi):
            os.mkdir(out_dir_snvfi)

        sub = samples_bulk.index(sample) + 1
        con = samples_bulk.index(bulk) + 1

        if sex is "M":
            ini_text = """SNV={0}\nSUB={1}\nCON={2}\nOUT_DIR={3}\nBLACKLIST=(\n'/hpc/pmc_vanboxtel/data/Mutation_blacklists/MSC_healthyBM_raw_variants.vcf'\n'/hpc/pmc_vanboxtel/data/dbSNP/nocosmic_dbsnp_137.b37.vcf'\n);\n\nSUB_GQ=99\nCON_GQ=10\nMQ=60\nQUAL=50\nCOV=20\nVAF=0.1\nFILTER=PASS\nCHR="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"\nCHR_NAM=NoY\nSNV_only=YES\nmax_alleles=2\nMAIL=F.M.Manders-2@prinsesmaximacentrum.nl\nCLEANUP=YES""".format(vcf, sub, con, out_dir_snvfi)
            ini_fname = out_dir + "SNVFI_" + sample + "_v1.2.ini"
            if not os.path.isfile(ini_fname):
                with open(ini_fname, "w") as ini:
                    ini.write(ini_text)

                command = "bash {0} {1} {2}".format(snvfi, snvfi_config, ini_fname)
                snvfi_bash_fname = out_dir_snvfi + "run_snvfi.sh"

                snvfi_bash = "#!/bin/bash\n\n#$ -S /bin/bash\n#$ -l h_rt=02:00:00\n#$ -cwd\n#$ -o {0}runsnvfi_output.txt\n#$ -e {0}runsnvfi_errors.txt\n#$ -M F.M.Manders-2@prinsesmaximacentrum.nl\n#$ -m as\n\nsource {1}\n\n{2}".format(out_dir_snvfi, guix_profile_source, command)
                with open(snvfi_bash_fname, "w") as snvfi_bash_file:
                    snvfi_bash_file.write(snvfi_bash)
                os.system("qsub {0}".format(snvfi_bash_fname))

                #os.system(command)

            #Filter on the X chromosome seperately with a high vaf, low depth and low sample gq.
            out_dir_snvfi = out_dir + sample + "_X/"
            if not os.path.isdir(out_dir_snvfi):
                os.mkdir(out_dir_snvfi)

            ini_text = """SNV={0}\nSUB={1}\nCON={2}\nOUT_DIR={3}\nBLACKLIST=(\n'/hpc/pmc_vanboxtel/data/Mutation_blacklists/MSC_healthyBM_raw_variants.vcf'\n'/hpc/pmc_vanboxtel/data/dbSNP/nocosmic_dbsnp_137.b37.vcf'\n);\n\nSUB_GQ=10\nCON_GQ=10\nMQ=60\nQUAL=50\nCOV=10\nVAF=0.99\nFILTER=PASS\nCHR="X"\nCHR_NAM=X\nSNV_only=YES\nmax_alleles=2\nMAIL=F.M.Manders-2@prinsesmaximacentrum.nl\nCLEANUP=YES""".format(vcf, sub, con, out_dir_snvfi)
            ini_fname = out_dir + "SNVFI_" + sample + "_X_v1.2.ini"
            if not os.path.isfile(ini_fname):
                with open(ini_fname, "w") as ini:
                    ini.write(ini_text)

                command = "bash {0} {1} {2}".format(snvfi, snvfi_config, ini_fname)
                snvfi_bash_fname = out_dir_snvfi + "run_snvfi.sh"

                snvfi_bash = "#!/bin/bash\n\n#$ -S /bin/bash\n#$ -l h_rt=02:00:00\n#$ -cwd\n#$ -o {0}runsnvfi_output.txt\n#$ -e {0}runsnvfi_errors.txt\n#$ -M F.M.Manders-2@prinsesmaximacentrum.nl\n#$ -m as\n\nsource {1}\n\n{2}".format(out_dir_snvfi, guix_profile_source, command)
                with open(snvfi_bash_fname, "w") as snvfi_bash_file:
                    snvfi_bash_file.write(snvfi_bash)
                os.system("qsub {0}".format(snvfi_bash_fname))

        elif sex is "F":
            ini_text = """SNV={0}\nSUB={1}\nCON={2}\nOUT_DIR={3}\nBLACKLIST=(\n'/hpc/pmc_vanboxtel/data/Mutation_blacklists/MSC_healthyBM_raw_variants.vcf'\n'/hpc/pmc_vanboxtel/data/dbSNP/nocosmic_dbsnp_137.b37.vcf'\n);\n\nSUB_GQ=99\nCON_GQ=10\nMQ=60\nQUAL=50\nCOV=20\nVAF=0.1\nFILTER=PASS\nCHR="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"\nCHR_NAM=NoY\nSNV_only=YES\nmax_alleles=2\nMAIL=F.M.Manders-2@prinsesmaximacentrum.nl\nCLEANUP=YES""".format(vcf, sub, con, out_dir_snvfi)
            ini_fname = out_dir + "SNVFI_" + sample + "_v1.2.ini"
            if not os.path.isfile(ini_fname):
                with open(ini_fname, "w") as ini:
                    ini.write(ini_text)

                command = "bash {0} {1} {2}".format(snvfi, snvfi_config, ini_fname)
                snvfi_bash_fname = out_dir_snvfi + "run_snvfi.sh"

                snvfi_bash = "#!/bin/bash\n\n#$ -S /bin/bash\n#$ -l h_rt=02:00:00\n#$ -cwd\n#$ -o {0}runsnvfi_output.txt\n#$ -e {0}runsnvfi_errors.txt\n#$ -M F.M.Manders-2@prinsesmaximacentrum.nl\n#$ -m as\n\nsource {1}\n\n{2}".format(out_dir_snvfi, guix_profile_source, command)
                with open(snvfi_bash_fname, "w") as snvfi_bash_file:
                    snvfi_bash_file.write(snvfi_bash)
                os.system("qsub {0}".format(snvfi_bash_fname))
        else:
            raise ValueError("Please provide a sex for fetus: {0} in the form of either M or F".format(fetus))



#Run indelfi per sample
        out_dir_indelfi = out_dir + sample + "_indel/"
        if not os.path.isdir(out_dir_indelfi):
            os.mkdir(out_dir_indelfi)


        indel_vcf = out_dir_indelfi + sample + "_" + bulk + "_QUAL250_RD10_VAF0.1_GQ10_flank100_MQ60_INDELs_autosomal_X_noEvidenceControl_noBlacklist.vcf"
        if not os.path.isfile(indel_vcf):
            sub_fromstart = sub + 8
            con_fromstart = con + 8
            command = "perl {0} -i {1} -s {2} -c {3} -o {4} -RD 10 -GQ 10 -bl /hpc/pmc_vanboxtel/data/Mutation_blacklists/MSC_healthyBM_raw_variants.vcf".format(indelfi, vcf, sub_fromstart, con_fromstart, out_dir_indelfi)
            indelfi_bash_fname = out_dir_indelfi + "run_indelfi.sh"
            indelfi_bash = "#!/bin/bash\n\n#$ -S /bin/bash\n#$ -l h_rt=00:30:00\n#$ -cwd\n#$ -o {0}runindelfi_output.txt\n#$ -e {0}runindelfi_errors.txt\n#$ -M F.M.Manders-2@prinsesmaximacentrum.nl\n#$ -m as\n\nsource {1}\n\n{2}".format(out_dir_indelfi, guix_profile_source, command)
            with open(indelfi_bash_fname, "w") as indelfi_bash_file:
                indelfi_bash_file.write(indelfi_bash)
            os.system("qsub {0}".format(indelfi_bash_fname))
            print("Ran infdelfi for sample {0} in fetus {1}".format(sample, fetus))

#Get callable regions of each sample and merge them with the callable regions of the bulk
        #GEt callable region
        bed_fname = dir_in + sample + "/" + sample + "_CallableLoci.bed"
        bed_callable_fname = out_callable + sample + "CallableLoci_CALLABLE.bed"
        if not os.path.isfile(bed_fname):
            bed_fname = dir_in + sample + "_dedup.realigned/" + sample + "_dedup.realigned_CallableLoci.bed"
        if os.path.isfile(bed_fname) and not os.path.isfile(bed_callable_fname):
            command = "grep 'CALLABLE' {0} > {1}".format(bed_fname, bed_callable_fname)
            os.system(command)
            print("Created CALLABLE file for sample {0} in fetus {1}".format(sample, fetus))

        #Intersect bed sample with bed bulk
        bed_intersected = out_callable + sample + "_" + bulk + "_CallableLoci.bed"
        if os.path.isfile(bed_callable_fname) and os.path.isfile(bulk_callable_fname) and not os.path.isfile(bed_intersected):
            command = "intersectBed -a {0} -b {1} > {2}".format(bed_callable_fname, bulk_callable_fname, bed_intersected)
            os.system(command)
            print("Intersected Callable bed file for sample {0} in fetus {1} with the bulk {2}".format(sample, fetus, bulk))

        #Sort the bedfile
        bed_sorted = out_callable + sample + "_" + bulk + "_CallableLoci_sorted.bed"
        if os.path.isfile(bed_intersected) and not os.path.isfile(bed_sorted):
            command = "sortBed -i {0} > {1}".format(bed_intersected, bed_sorted)
            os.system(command)
            print("Sorted intersected bed for sample {0} in fetus {1}".format(sample, fetus))

        #Merge overlapping regions of the bedfile
        bed_merged = out_callable + sample + "_" + bulk + "_CallableLoci_merged.bed"
        if os.path.isfile(bed_sorted) and not os.path.isfile(bed_merged):
            command = "mergeBed -i {0} > {1}".format(bed_sorted, bed_merged)
            os.system(command)
            print("Merged overlapping regions in bed for sample {0} in fetus {1}".format(sample, fetus))

        #Filter for autosomal regions
        bed_merged_autosomal = out_callable + sample + "_" + bulk + "_CallableLoci_merged_autosomal.bed"
        if os.path.isfile(bed_merged) and not os.path.isfile(bed_merged_autosomal):
            command = "sed '/[XYMT]/d' {0} > {1}".format(bed_merged, bed_merged_autosomal)
            os.system(command)
            print("Filtered for autosomal sites for sample {0} in fetus {1}".format(sample, fetus))

        #Count surveyed region of the genome
        surveyed_region = out_callable + sample + "_" + bulk + "_surveyed.txt"
        if os.path.isfile(bed_merged_autosomal) and not os.path.isfile(surveyed_region):
            command = "cat " + bed_merged_autosomal + " | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' > " + surveyed_region
            os.system(command)
            print("Counted surveyed region (autosomal) for sample {0} in fetus {1}".format(sample, fetus))

        #Filter for autosomal + X chromosome regions
        bed_merged_autosomalX = out_callable + sample + "_" + bulk + "_CallableLoci_merged_autosomal_X.bed"
        if os.path.isfile(bed_merged) and not os.path.isfile(bed_merged_autosomalX):
            command = "sed '/[YMT]/d' {0} > {1}".format(bed_merged, bed_merged_autosomalX)
            os.system(command)
            print("Filtered for autosomal + X sites for sample {0} in fetus {1}".format(sample, fetus))

        #Count surveyed region of the genome
        surveyed_regionX = out_callable + sample + "_" + bulk + "_surveyed_autosomal_X.txt"
        if os.path.isfile(bed_merged_autosomalX) and not os.path.isfile(surveyed_regionX):
            command = "cat " + bed_merged_autosomalX + " | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' > " + surveyed_regionX
            os.system(command)
            print("Counted surveyed region (autosomal+X) for sample {0} in fetus {1}".format(sample, fetus))
