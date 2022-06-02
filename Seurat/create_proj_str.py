import os
import sys
import argparse


'''Usage:create_project_str.py [-h] -p <project path>'''
'''create a list of the directory names that need to be created for a project'''
dir_list = ['Quality_figures', 'Seurat_objects', 'Integration_results', 'DE_Analysis', 'Enrichment_files','data']

def msg():
    return ''' This script checks if a project directory and sub-directory exists,
        and creates them if necessary. Please provide a project location. Exiting.... '''


def parse_arguments():
    '''Adding the command line arguments to be provided while running the script'''
    parser = argparse.ArgumentParser(description='Process command line arguments.', usage=msg())
    parser.add_argument('-p', '--path',help ='project path')
    args = parser.parse_args()
    return args


def check_proj_path():

    '''This function checks if main project directory exists...'''
    input_path = parse_arguments()
    dir_path = input_path.path
    if os.path.isdir(dir_path):
        #print(f"directory {dir_path} exists!")
        return dir_path
    else:
        print("Project dir: " + dir_path + " does not exist")
        sys.exit()


def copy_config(new_path):

    '''Copy the config file to the user-provided sub-directory'''
    cwd = os.getcwd()
    config_file = "config.yml"
    config_path = os.path.join(cwd, config_file)
    print("Config file: " + config_path)

    if os.path.isfile(config_path):
        try:
            cmd = "cp " +  config_path  + " " + new_path
            print(cmd)
            return_val = os.system(cmd)
            #print("Return value: " + str(return_val))
            if (return_val == 0):
                print("Finished copying the config file to the " + new_path)
        except:
            print("Return value: " + str(return_val))
            print("Failed to copy the configuration file to the " +  new_path)


def create_dir(new_path):
    '''A function to check if the directory exists & create one if necessary & copy the config file'''

    if os.path.isdir(new_path):
        print("Sub-directory, {new_path} already exists!")
    else:
        try:
            os.mkdir(new_path)
            print("Finished creating the" + new_path)
        except:
            print("Failed to create the subdirectory " + new_path)
        else:
	        pass

################
# Main function:  This function handles menu and flow control.
################
def main():

    proj_dir = check_proj_path()
    print("Project Location: " + proj_dir)

    for sub_directory in dir_list:
        new_path = os.path.join(proj_dir, sub_directory)
        print("\nChecking/Creating sub-directory: " + new_path)
        create_dir(new_path)

        if(sub_directory == "Spreadsheets"):
            print("\nCopying the config.yml file")
            copy_config(new_path)
        else:
            pass

if __name__ == '__main__':
    main()
