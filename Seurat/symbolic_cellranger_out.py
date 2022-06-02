import os
import sys
import argparse
import logging
import pandas as pd
import shutil

'''Usage:create_cellranger_output.py -cp <path containing cellranger output for all samples> -m csv file containing sample name and condition information'''

def msg():
    return ''' This script create the directory structure need to run Seurat R script. Exiting.... '''

def parse_arguments():
    '''Adding the command line arguments to be provided while running the script'''
    parser = argparse.ArgumentParser(description='Process command line arguments.', usage=msg())
    parser.add_argument('-cp', '--cellranger_path', required=True, help='cellranger output path')
    parser.add_argument('-m', '--metadata', metavar='FILE', required=True, help='metadata file')
    args = parser.parse_args()
    return args

cmd_line_args = parse_arguments()
cellranger_path = cmd_line_args.cellranger_path
metadata = cmd_line_args.metadata
print("User provided folder path: " + cellranger_path + " and metadata info: " + metadata )


def make_sample_directories():


    #Iterate through the sample list
    if os.path.isdir(cellranger_path):
        print("Received path to cell ranger output: " + cellranger_path)

    if os.path.isfile(metadata):
        print("Received metadata file: " + metadata)

    metadata_file = pd.read_csv(metadata)
    sample_names=metadata_file['orig.ident']

    for sample in sample_names:
        os.mkdir('data/' + sample)
        target_dir='data/' + sample + "/"
        source_dir = cellranger_path + sample + '/outs/filtered_feature_bc_matrix/'
    
    for filename in os.listdir(source_dir):
        shutil.copy(source_dir + filename, target_dir)
    
make_sample_directories()
