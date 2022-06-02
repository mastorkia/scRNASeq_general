import os
import sys
import argparse
import logging

'''Usage:run_cellranger.py -f <file list containing raw fastq file(s) location CIM.list> -c <config file> -a <purpose>'''
#python run_cellranger.py -f file.list -c config_cellranger.yml -a 'MULTI'
#cat file.list
#/share/data/RNA_Seq/10X/Raw_fastq/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_I1_001.fastq.gz
#/share/data/RNA_Seq/10X/Raw_fastq/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz
#/share/data/RNA_Seq/10X/Raw_fastq/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz

def msg():
    return ''' This script runs cellranger program from totsalseq platform to map the single cell reads, demultiplex and generate the count
                matrices. Exiting.... '''

def parse_arguments():
    '''Adding the command line arguments to be provided while running the script'''
    parser = argparse.ArgumentParser(description='Process command line arguments.', usage=msg())
    parser.add_argument('-f', '--file_list', metavar='FILE', required=True, help='fastq file list')
    parser.add_argument('-c', '--config', metavar='FILE', required=True, help='config file')
    parser.add_argument('-a', '--aim', metavar='MULTI/COUNT', required=True, help='MULTI/COUNT')
    args = parser.parse_args()
    return args

cmd_line_args = parse_arguments()
file_list = cmd_line_args.file_list
print("User provided file list: " + file_list)
aim = cmd_line_args.aim
print("The aim of the current analysis is: " + aim)


def make_config_hash():

    print("Reading the config file")

    config_dict = {}
    config_file = cmd_line_args.config
    if os.path.isfile(config_file):
        print("Received config file: " + config_file)
    else:
        sys.exit()

    config = open( config_file, 'r')

    for line in config:
        this_line = line.strip()
        this_key = this_line.split(":")[0]
        this_val = this_line.split(":")[1]
        config_dict[this_key] = this_val

    return config_dict


def make_sample_dict():

    sample_dict = {}

    #Iterate through the file list
    if os.path.isfile(file_list):
        print("Received fastq file list: " + file_list)

    fastq_list = open(file_list, 'r')

    for fastq in fastq_list:
        this_fastq = fastq.strip()
        if os.path.isfile(this_fastq):
            print("Found: " + this_fastq)
        else:
            print("this_fastq does not exist. Please check the file")

        dir_name = os.path.dirname(this_fastq)
        file_name = os.path.basename(this_fastq)
        sample_name = file_name.split("_")[0]
        sample_dict[sample_name] = dir_name

    return sample_dict


def run_cellranger():

    my_config_dict = make_config_hash()
    print("Printing config hash")
    print(my_config_dict)

    sample_dictionary = make_sample_dict()
    print("Printing sample dictionary")
    print(sample_dictionary)

    for sample_name, dir_name in sample_dictionary.items():
        this_transcriptome_loc = my_config_dict['transcriptome_loc']
        print("Transcriptomic Reference Location: "  + this_transcriptome_loc)
        cellranger_loc = my_config_dict['software_loc']
        print("Software location: " + cellranger_loc)
        expected_cells = my_config_dict['expect_cells']
        num_cores = my_config_dict['local_cores']
        mem = my_config_dict['local_mem']
        print("Expected Cells: " + expected_cells)
        print("Number of Cores: " + num_cores)
        print("Memory Assigned: " + mem)
        job_name = sample_name + "_cellranger"
        print(job_name)

        if aim == 'MULTI':
            cmd = "echo " + cellranger_loc + " multi --id=" + sample_name + " --csv=" + sample_name +".config" \
            " | qsub -S /bin/sh -V -N " + job_name + " -l h_vmem=" + mem + " -pe smp " +  num_cores + " -e " + dir_name + "/" + " -o " +  dir_name + "/" +  " -j y -cwd"

        elif aim == 'COUNT':
            cmd = "echo " + cellranger_loc + " count --id=" + sample_name + " --libraries=" + sample_name + "_count_Abs.config " + "--feature-ref=" + dir_name + "/CMO_reference.csv --expect-cells=" + expected_cells + " --transcriptome=" + this_transcriptome_loc + \
            " | qsub -S /bin/sh -V -N " + job_name + " -l h_vmem=" + mem + " -pe smp " + num_cores + " -e " + dir_name + "/" + " -o " + dir_name + "/" + " -j y -cwd"

        #generate a SLURM command
        print(cmd)

        os.chdir(dir_name)
        current_dir = os.getcwd()
        print("\nWorking dir: " + current_dir)
        return_val = os.system(cmd)

        if return_val == 0:
            print("Cellranger job successfully submitted for sample:  " + sample_name + " and purpose " + aim)

        else:
            print("Return value: " + str(return_val))
            print("Cellranger job cannot be successfully submitted for sample: " + sample_name)

run_cellranger()
