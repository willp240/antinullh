import os
import argparse
import configparser
import json

def check_dir(dname):
    """Check if directory exists, create it if it doesn't"""
    if(dname[-1] != "/"):
        dname = dname + "/"
    direc = os.path.dirname(dname)
    try:
        os.stat(direc)
    except:
        os.makedirs(direc)
        print ("Made directory %s...." % dname)
    return dname


def pycondor_submit(job_name, exec_name, out_dir, run_dir, env_file, fit_config, rates_config, pdf_config, syst_config, osc_config, walltime, sleep_time = 1, priority = 5):
    '''
    submit a job to condor, write a sh file to source environment and execute command
    then write a submit file to be run by condor_submit
    '''

    print (job_name)
    if exec_name == "make_osc_grids":
        batch_name, index = job_name.split('_index_')
    else:
        index = ""
        batch_name = job_name

    ### set a condor path to be called later
    
    condor_path = "{0}/".format(out_dir)
    exec_path = run_dir + "/bin/" + exec_name

    configs_path = os.path.abspath('{0}/cfg'.format(condor_path))
    check_dir(configs_path)

    if fit_config != "":
        os.system("cp " + str(fit_config) + " " + str(configs_path) + "/" + os.path.basename(fit_config) )
        fit_config = configs_path + "/" + fit_config.split("/")[-1]
    if rates_config != "":
        os.system("cp " + str(rates_config) + " " + str(configs_path) + "/" + os.path.basename(rates_config) )
        rates_config = configs_path + "/" + rates_config.split("/")[-1]
    if pdf_config != "":
        os.system("cp " + str(pdf_config) + " " + str(configs_path) + "/" + os.path.basename(pdf_config) )
        pdf_config = configs_path + "/" + pdf_config.split("/")[-1]
    if syst_config != "":
        os.system("cp " + str(syst_config) + " " + str(configs_path) + "/" + os.path.basename(syst_config) )
        syst_config = configs_path + "/" + syst_config.split("/")[-1]
    if osc_config != "":
        os.system("cp " + str(osc_config) + " " + str(configs_path) + "/" + os.path.basename(osc_config) )
        osc_config = configs_path + "/" + osc_config.split("/")[-1]

    other_commands = 'sleep $[($RANDOM%' + str(sleep_time+1) + ')+1]s'

    out_macro_text = "#!/usr/bin/sh  \n" + \
                     "source " + env_file + "\n" + \
                     "source " + run_dir + "/env.sh " + "\n" + \
                     "echo $HOME \n" + \
                     "cd " + str(run_dir) + "\n" + \
                     str(other_commands) + "\n" + \
                     str(exec_path) + " " + str(fit_config) + " " + str(rates_config) + " " + str(pdf_config) + " " + str(syst_config) + " " + str(osc_config) + " " + index

    sh_filepath = "{0}/sh/".format(condor_path) + str(job_name).replace("/", "") + '.sh'
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
    getenv = "True"

    submit_filepath = os.path.join(submit_path, job_name)
    submit_filepath += ".submit"
    out_submit_text = "executable              = " + str(sh_filepath) + "\n" + \
                     "universe                 = " + str(universe) + "\n" + \
                     "output                   = " + str(output_path) + "/" + str(job_name) + ".output\n" + \
                     "error                    = " + str(error_path) + "/" + str(job_name) + ".error\n" + \
                     "log                      = " + str(log_path) + "/" + str(job_name) + ".log\n" + \
                     "notification             = " + str(notification) + "\n" + \
                     "priority                 = " + str(priority) + "\n" + \
                     "getenv                   = " + str(getenv) + "\n" + \
                     "allowed_execute_duration = " + str(walltime) + " \n" + \
                     "queue "+str(n_rep)+"\n"

    ## check and create output path
    if not os.path.exists(os.path.dirname(submit_filepath)):
        os.makedirs(os.path.dirname(submit_filepath))
    out_submit_file = open(submit_filepath, "w")
    out_submit_file.write(out_submit_text)
    out_submit_file.close()


    command = 'condor_submit -batch-name \"' + batch_name +'\" ' + submit_filepath
    print ("executing job: " + command)
    os.system(command)


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Submit jobs to condor")
    parser.add_argument('exe', type=str, help='executable')
    parser.add_argument('out_dir', type=str, help='output directory')
    parser.add_argument('-r', "--run_dir", type=str, default=".", help='directory to run executable from')
    parser.add_argument('-e', "--env_file", type=str, default="", help='environment variable file to source')
    parser.add_argument('-f', "--fit_cfg", type=str, default="", help='fit config path')
    parser.add_argument('-i', "--rates_cfg", type=str, default="", help='event rates config path')
    parser.add_argument('-p', "--pdf_cfg", type=str, default="", help='pdf config path')
    parser.add_argument('-s', "--syst_cfg", type=str, default="", help='syst config path')
    parser.add_argument('-o', "--osc_cfg", type=str, default="", help='osc grid config path')
    parser.add_argument("-n", "--num_jobs", type=int, default=1, help="how many identical jobs would you like to run?")
    parser.add_argument("-w", "--wall_time", type=int, default=86400, help="what's the maximum runtime (in seconds, default 1 day)?")
    args = parser.parse_args()

    ## check if output and condor directories exist, create if they don't
    exec_name = args.exe
    out_dir = check_dir(args.out_dir)
    base_name = out_dir.split("/")[-2] #os.path.basename(out_dir)
    run_dir = args.run_dir
    env_file = args.env_file
    fit_config = args.fit_cfg
    rates_config = args.rates_cfg
    pdf_config = args.pdf_cfg
    syst_config = args.syst_cfg
    osc_config = args.osc_cfg
    walltime = args.wall_time
    if fit_config != "":
        fit_config = run_dir + "/" + args.fit_cfg
    if rates_config != "":
        rates_config = run_dir + "/" + args.rates_cfg
    if pdf_config != "":
        pdf_config = run_dir + "/" + args.pdf_cfg
    if syst_config != "":
        syst_config = run_dir + "/" + args.syst_cfg
    if osc_config != "":
        osc_config = run_dir + "/" + args.osc_cfg

    if exec_name == "make_osc_grids":
        # Get reactors json from oscgrid
        config = configparser.ConfigParser()
        config.read(osc_config)

        # Access the reactors json filepath under the 'summary' section
        reactorsjson_path = config.get('summary', 'reactorsjson')

        # Load the JSON file
        with open(reactorsjson_path) as f:
            data = json.load(f)

        # Loop over each key-value pair
        for index, values in data.items():
            # Convert the key to an integer
            int_index = int(index)
            # Get the string and double values from the list
            reac_name, distance_val = values

            job_name = base_name + "_index_{0}".format(int_index)

            log_dir = check_dir("{0}/log/".format(out_dir))
            error_dir = check_dir("{0}/error/".format(out_dir))
            sh_dir = check_dir("{0}/sh/".format(out_dir))
            submit_dir = check_dir("{0}/submit/".format(out_dir))
            output_dir = check_dir("{0}/output/".format(out_dir))

            pycondor_submit(job_name, exec_name, out_dir, run_dir, env_file, fit_config, rates_config, pdf_config, syst_config, osc_config, walltime, sleep_time = 1, priority = 5)

    else:
        for i in range(args.num_jobs):

            job_name = base_name + "_{0}".format(i)

            log_dir = check_dir("{0}/log/".format(out_dir))
            error_dir = check_dir("{0}/error/".format(out_dir))
            sh_dir = check_dir("{0}/sh/".format(out_dir))
            submit_dir = check_dir("{0}/submit/".format(out_dir))
            output_dir = check_dir("{0}/output/".format(out_dir))

            pycondor_submit(job_name, exec_name, out_dir, run_dir, env_file, fit_config, rates_config, pdf_config, syst_config, osc_config, walltime, sleep_time = 1, priority = 5)
