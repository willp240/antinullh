import sys
import os
import string
import argparse

def check_dir(dname):
    """Check if directory exists, create it if it doesn't"""
    if(dname[-1] != "/"):
        dname = dname + "/"
    direc = os.path.dirname(dname)
    try:
        os.stat(direc)
    except:
        os.makedirs(direc)
        print "Made directory %s...." % dname
    return dname

def pycondor_submit(job_batch, job_id, run_path, dataset_dir, dataset_name, fit_dir, dim, fit_config, cut_config, pdf_config, syst_config, env_file, sleep_time = 1, priority = 5):
    '''
    submit a job to condor, write a sh file to source environment and execute command
    then write a submit file to be run by condor_submit
    '''

    print job_id

    ### set a condor path to be called later
    
    condor_path = "{0}/".format(fit_dir)
    exec_path = run_path + "/bin/fit_dataset"

    configs_path = os.path.abspath('{0}/cfg'.format(condor_path))
    check_dir(configs_path)
    os.system("cp " + str(fit_config) + " " + str(configs_path) + "/" + os.path.basename(fit_config) )
    os.system("cp " + str(pdf_config) + " " + str(configs_path) + "/" + os.path.basename(pdf_config) )
    os.system("cp " + str(cut_config) + " " + str(configs_path) + "/" + os.path.basename(cut_config) )
    os.system("cp " + str(syst_config) + " " + str(configs_path) + "/" + os.path.basename(syst_config) )


    dataset_file = dataset_dir + dataset_name

    other_commands = 'sleep $[($RANDOM%' + str(sleep_time+1) + ')+1]s'

    out_macro_text = "#!/usr/bin/sh  \n" + \
                     "source " + str(env_file) + "\n" + \
                     "cd " + str(run_path) + "\n" + \
                     str(other_commands) + "\n" + \
                     str(exec_path) + " " + str(configs_path) + "/" + os.path.basename(fit_config) + " " + str(configs_path) + "/" + os.path.basename(pdf_config) + " " + str(configs_path) + "/" + os.path.basename(cut_config) + " " + str(configs_path) + "/" + os.path.basename(syst_config) + " " + str(dataset_file) + " " + str(dim) + " " + str(fit_dir) + "\n"

    sh_filepath = "{0}/sh/".format(condor_path) + str(job_id).replace("/", "") + '.sh'
    if not os.path.exists(os.path.dirname(sh_filepath)):
        os.makedirs(os.path.dirname(sh_filepath))
    sh_file = open(sh_filepath, "w")
    sh_file.write(out_macro_text)
    sh_file.close()
    os.chmod(sh_filepath, 0o777)
    
    ### now create submit file
    error_path = os.path.abspath('{0}/error'.format(condor_path))
    output_path = os.path.abspath('{0}/output'.format(condor_path))
    log_path = os.path.abspath('{0}/log'.format(condor_path))
    submit_path = os.path.abspath('{0}/submit'.format(condor_path))

    universe = "vanilla"
    notification = "never"
    n_rep = 1
    getenv = "False"

    submit_filepath = os.path.join(submit_path, job_id)
    submit_filepath += ".submit"
    out_submit_text = "executable              = " + str(sh_filepath) + "\n" + \
                     "universe                = " + str(universe) + "\n" + \
                     "output                  = " + str(output_path) + "/" + str(job_id) + ".output\n" + \
                     "error                   = " + str(error_path) + "/" + str(job_id) + ".error\n" + \
                     "log                     = " + str(log_path) + "/" + str(job_id) + ".log\n" + \
                     "notification            = " + str(notification) + "\n" + \
                     "priority                = " + str(priority) + "\n" + \
                     "getenv                  = " + str(getenv) + "\n" + \
                     "queue "+str(n_rep)+"\n"

    ## check and create output path
    if not os.path.exists(os.path.dirname(submit_filepath)):
        os.makedirs(os.path.dirname(submit_filepath))
    out_submit_file = open(submit_filepath, "w")
    out_submit_file.write(out_submit_text)
    out_submit_file.close()

    command = 'condor_submit -batch-name \"' + job_batch+'\" ' + submit_filepath
    print "executing job: " + command
    os.system(command)


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Launch a load of fits")
    parser.add_argument('dataset_dir', type=str, help='Dataset directory. Top level directory where ntuples e_r_class_ntups directories have been written in. Intended structure is to write all youur datasets made from these ntuples into this directory, as well as any fits to those datasets. ')
    parser.add_argument('dataset_name', type=str, help='Dataset filename. Should be an h5 file within dataset_dir (with extension)')
    parser.add_argument('fit_dir', type=str, help='Fit directory name. This will be made within dataset directory.')
    parser.add_argument('-d', "--dims", type=str, help='How any dimensions?')
    parser.add_argument('-f', "--fit_cfg", type=str, help='Fit config path')
    parser.add_argument('-c', "--cut_cfg", type=str, help='Cut config path')
    parser.add_argument('-p', "--pdf_cfg", type=str, help='PDF config path')
    parser.add_argument('-s', "--syst_cfg", type=str, help='Systematic config path')
    parser.add_argument('-e', "--env_file", type=str, help='Environment file to source. Absolute path.')
    parser.add_argument("-n", "--no_sims", type=int,
                       default=10,
                       help="How many identical fits would you like to run?")
    args = parser.parse_args()

    ## check if output and condor directories exist, create if they don't
    bb_dir = os.getenv('BB_DIR')
    data_dir = os.getenv('DATA_DIR')
    dataset_dir = data_dir + args.dataset_dir
    dataset_dir = check_dir(dataset_dir)
    dataset_name = args.dataset_name
    base_name = args.fit_dir
    over_dir = check_dir(dataset_dir+"/"+base_name) 
    fit_config = bb_dir + "../results/" + args.fit_cfg
    cut_config = bb_dir + "../cuts/" + args.cut_cfg
    pdf_config = bb_dir + "../pdfs/" + args.pdf_cfg
    syst_config = bb_dir + "../systs/" + args.syst_cfg
    env_file = args.env_file

    dims = (args.dims)

    for i in range(args.no_sims):

        job_name = base_name + "_{0}".format(i)
        fit_dir = over_dir + "/" + job_name
        fit_dir = check_dir(fit_dir)
        dims = (args.dims)
        
        condor_directory = "{0}/condor".format(bb_dir)
        
        log_dir = check_dir("{0}/log/".format(fit_dir))
        error_dir = check_dir("{0}/error/".format(fit_dir))
        sh_dir = check_dir("{0}/sh/".format(fit_dir))
        submit_dir = check_dir("{0}/submit/".format(fit_dir))
        output_dir = check_dir("{0}/output/".format(fit_dir))

        job_id = "{0}_{1}".format(base_name,i)
        batch_id = "{0}".format(args.fit_dir)
        pycondor_submit(batch_id, job_id, bb_dir, dataset_dir, dataset_name, fit_dir, dims, fit_config, cut_config, pdf_config, syst_config, env_file, sleep_time = 1, priority = 5)
