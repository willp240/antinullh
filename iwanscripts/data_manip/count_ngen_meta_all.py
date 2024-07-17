'''
Script to count ngen for all event types in a production.
'''

import os
import subprocess
import string 

ENVIRONMENT_PATH = "/home/kroupova/env_git_rat.sh"
SUBMIT_STRING = "qsub -l cput=1:59:59"
COUNT_SCRIPT = "/home/kroupova/bb_sigex/data_manip/count_ngen_meta.py"


def write_shell_script(subdirectory, name, outdir, shell_script_name):
    
    with open(shell_script_name, "w") as f:
        f.write("""source {0}
cd {1}
python {2} {2} >> counts_{3}.txt
        """.format(ENVIRONMENT_PATH, COUNT_SCRIPT, outdir, subdirectory+"/*", name))
    
def submit_job(shell_script_name):
    subprocess.check_call("{0} {1}".format(SUBMIT_STRING, shell_script_name), shell = True)
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("production_directory", type = str)
    parser.add_argument("outdir", type=str)
    args = parser.parse_args()

    for subdir in os.listdir(args.production_directory):
        print args.production_directory+subdir
        shell_name = subdir+".sh"
        write_shell_script(os.path.join(args.production_directory, subdir), subdir, args.outdir, shell_name)
        submit_job(shell_name)
