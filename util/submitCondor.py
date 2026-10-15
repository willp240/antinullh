import os
import argparse
import re
import shlex
import shutil
import subprocess

def check_dir(dname):
    """Check if directory exists, create it if it doesn't"""
    os.makedirs(dname, exist_ok=True)
    return os.path.join(dname, "")


def copy_config(path, configs_path):
    if not path:
        return ""
    destination = os.path.join(configs_path, os.path.basename(path))
    shutil.copy2(path, destination)
    return destination


def update_output_directory(fit_config, out_dir):
    """Update the copied config's summary field, preserving other config text."""
    with open(fit_config, "r") as file:
        lines = file.readlines()
    section = None
    matches = []
    for index, line in enumerate(lines):
        heading = re.match(r"\s*\[([^]]+)\]", line)
        if heading:
            section = heading.group(1).strip()
        elif section == "summary" and re.match(r"\s*output_directory\s*=", line):
            matches.append(index)
    if len(matches) != 1:
        raise ValueError(f"{fit_config}: expected exactly one output_directory in [summary]")
    lines[matches[0]] = f"output_directory = {out_dir}\n"
    with open(fit_config, "w") as file:
        file.writelines(lines)


def pycondor_submit(job_name, exec_name, out_dir, run_dir, env_file, fit_config, event_config, pdf_config, syst_config, walltime, mem, sleep_time = 1, priority = 5, batch_name = ""):
    '''
    submit a job to condor, write a sh file to source environment and execute command
    then write a submit file to be run by condor_submit
    '''

    print (job_name)
    if batch_name == "":
        batch_name = job_name

    out_dir = os.path.abspath(out_dir)
    run_dir = os.path.abspath(run_dir)
    if env_file:
        env_file = os.path.abspath(os.path.expanduser(env_file))
        if not os.path.isfile(env_file):
            raise FileNotFoundError(env_file)

    # Set a condor path to be called later
    condor_path = "{0}/".format(out_dir)
    exec_path = run_dir + "/bin/" + exec_name

    configs_path = os.path.abspath('{0}/cfg'.format(condor_path))
    check_dir(configs_path)

    fit_config = copy_config(fit_config, configs_path)
    event_config = copy_config(event_config, configs_path)
    pdf_config = copy_config(pdf_config, configs_path)
    syst_config = copy_config(syst_config, configs_path)

    # Update output_directory in the fit config, if there is one
    if fit_config != "":
        update_output_directory(fit_config, out_dir)

    if exec_name == "makeFixedOscTree":
        exec_args = [fit_config]
    elif exec_name == "prune_trees":
        exec_args = [event_config, fit_config]
    else:
        exec_args = [fit_config, event_config, pdf_config, syst_config]
    if not all(exec_args):
        raise ValueError(f"Missing required config argument for {exec_name}")

    # Write sh file to be ran
    commands = ["#!/bin/bash", "set -e"]
    if env_file:
        commands.append("source " + shlex.quote(env_file))
    commands.append("cd " + shlex.quote(run_dir))
    commands.append(f"sleep $((RANDOM % {sleep_time + 1} + 1))s")
    commands.append("exec " + shlex.join([exec_path] + exec_args))
    out_macro_text = "\n".join(commands) + "\n"

    sh_filepath = "{0}/sh/".format(condor_path) + str(job_name).replace("/", "") + '.sh'
    if not os.path.exists(os.path.dirname(sh_filepath)):
        os.makedirs(os.path.dirname(sh_filepath))
    sh_file = open(sh_filepath, "w")
    sh_file.write(out_macro_text)
    sh_file.close()
    os.chmod(sh_filepath, 0o755)
    
    # Now create submit file
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
    out_submit_text = "executable              = " + sh_filepath + "\n" + \
                     "universe                 = " + str(universe) + "\n" + \
                     "output                   = " + os.path.join(output_path, job_name + ".output") + "\n" + \
                     "error                    = " + os.path.join(error_path, job_name + ".error") + "\n" + \
                     "log                      = " + os.path.join(log_path, job_name + ".log") + "\n" + \
                     "notification             = " + str(notification) + "\n" + \
                     "priority                 = " + str(priority) + "\n" + \
                     "getenv                   = " + str(getenv) + "\n" + \
                     "allowed_execute_duration = " + str(walltime) + " \n" + \
                     "request_memory           = " + str(mem) + " \n" + \
                     "queue "+str(n_rep)+"\n"

    # Check and create output path
    if not os.path.exists(os.path.dirname(submit_filepath)):
        os.makedirs(os.path.dirname(submit_filepath))
    out_submit_file = open(submit_filepath, "w")
    out_submit_file.write(out_submit_text)
    out_submit_file.close()

    # Lez do dis
    command = ["condor_submit", "-batch-name", batch_name, submit_filepath]
    print ("executing job: " + shlex.join(command))
    subprocess.run(command, check=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Submit jobs to condor")
    parser.add_argument('exe', type=str, help='executable')
    parser.add_argument('out_dir', type=str, help='output directory')
    parser.add_argument('-r', "--run_dir", type=str, default=".", help='directory to run executable from')
    parser.add_argument('-e', "--env_file", type=str, default="", help='complete environment setup script to source; otherwise inherit the submission environment')
    parser.add_argument('-f', "--fit_cfg", type=str, default="", help='fit config path')
    parser.add_argument('-i', "--event_cfg", type=str, default="", help='event config path')
    parser.add_argument('-p', "--pdf_cfg", type=str, default="", help='pdf config path')
    parser.add_argument('-s', "--syst_cfg", type=str, default="", help='syst config path')
    parser.add_argument("-n", "--num_jobs", type=int, default=1, help="how many identical jobs would you like to run?")
    parser.add_argument("-w", "--wall_time", type=int, default=86400, help="what's the maximum runtime (in seconds, default 1 day)?")
    parser.add_argument("-m", "--mem", type=float, default=300, help="what's the maximum memory (in MB, default 300 MB)?")
    parser.add_argument("-j", "--job_name", type=str, default="", help='job name')
    args = parser.parse_args()

    # Check if output and condor directories exist, create if they don't
    exec_name = args.exe
    out_dir = check_dir(args.out_dir)
    base_name = out_dir.split("/")[-2]
    run_dir = os.path.abspath(args.run_dir)
    env_file = args.env_file
    fit_config = args.fit_cfg
    event_config = args.event_cfg
    pdf_config = args.pdf_cfg
    syst_config = args.syst_cfg
    walltime = args.wall_time
    mem = args.mem
    job_name = args.job_name
    if fit_config != "":
        fit_config = os.path.join(run_dir, args.fit_cfg)
    if event_config != "":
        event_config = os.path.join(run_dir, args.event_cfg)
    if pdf_config != "":
        pdf_config = os.path.join(run_dir, args.pdf_cfg)
    if syst_config != "":
        syst_config = os.path.join(run_dir, args.syst_cfg)

    if job_name == "":
        job_name = base_name

    # Submit independent jobs
    for i in range(args.num_jobs):

        this_job_name = job_name + "_{0}".format(i)

        # With several jobs, give each its own subdirectory, otherwise they share
        # one fit config output_directory and overwrite each other's results
        if args.num_jobs > 1:
            this_out_dir = check_dir("{0}/{1}/".format(out_dir, this_job_name))
        else:
            this_out_dir = out_dir

        log_dir = check_dir("{0}/log/".format(this_out_dir))
        error_dir = check_dir("{0}/error/".format(this_out_dir))
        sh_dir = check_dir("{0}/sh/".format(this_out_dir))
        submit_dir = check_dir("{0}/submit/".format(this_out_dir))
        output_dir = check_dir("{0}/output/".format(this_out_dir))

        pycondor_submit(this_job_name, exec_name, this_out_dir, run_dir, env_file, fit_config, event_config, pdf_config, syst_config, walltime, mem, sleep_time = 1, priority = 5, batch_name = job_name)
