import os
import argparse

numvalsu = 99
numvalsth = 99

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

def read_fitcfg(fit_config):
    th_min = None
    th_max = None
    u_min = None
    u_max = None

    with open(fit_config, "r") as file:
        lines = file.readlines()

    for iline, line in enumerate(lines):
        if line.startswith("[geonu_U_norm]"):
            for jline in range(iline + 1, len(lines)):
                if lines[jline].startswith("min ="):
                    u_min = lines[jline].split("=")[1].strip()
                    break
            for jline in range(iline + 1, len(lines)):
                if lines[jline].startswith("max ="):
                    u_max = lines[jline].split("=")[1].strip()
                    break

    for iline, line in enumerate(lines):
        if line.startswith("[geonu_Th_norm]"):
            for jline in range(iline + 1, len(lines)):
                if lines[jline].startswith("min ="):
                    th_min = lines[jline].split("=")[1].strip()
                    break
            for jline in range(iline + 1, len(lines)):
                if lines[jline].startswith("max ="):
                    th_max = lines[jline].split("=")[1].strip()
                    break

    return u_min, u_max, th_min, th_max

def pycondor_submit(job_name, exec_name, out_dir, run_dir, env_file,
                    fit_config, event_config, pdf_config, syst_config, osc_config,
                    walltime, mem, th, start_idx, end_idx,
                    sleep_time=1, priority=5):
    '''
    Submit a job to condor, write a sh file to source environment and execute command
    then write a submit file to be run by condor_submit
    '''

    print (job_name)
    batch_name = job_name.split('_', 1)[0]
    first_u = float(u_min) + start_idx*(float(u_max)-float(u_min))/numvalsu
    first_u = "{:.3f}".format(first_u)
    last_u = float(u_min) + (end_idx-1)*(float(u_max)-float(u_min))/numvalsu
    last_u = "{:.3f}".format(last_u)
    out_dir_base = out_dir.split("/")[-2] + "_Th{0}_U{1}-{2}".format(th, first_u, last_u)

    condor_path = "{0}/".format(out_dir)
    exec_path = run_dir + "/bin/" + exec_name

    configs_path = os.path.abspath('{0}/th{1}/cfg'.format(condor_path, th))
    check_dir(configs_path)

    with open(fit_config, "r") as file:
        lines = file.readlines()

    # Loop over U values, write fit_config with U and Th values in
    for iFit in range(start_idx, end_idx):
        u = float(u_min) + iFit*(float(u_max)-float(u_min))/numvalsu
        u = "{:.3f}".format(u)

        # Process the file and make updates
        for iline, line in enumerate(lines):
            if line.startswith("[geonu_Th_norm]"):
                for jline in range(iline + 1, len(lines)):
                    if lines[jline].startswith("nom ="):
                        lines[jline] = f"nom = {th}\n"
                        break
            elif line.startswith("[geonu_U_norm]"):
                for jline in range(iline + 1, len(lines)):
                    if lines[jline].startswith("nom ="):
                        lines[jline] = f"nom = {u}\n"
                        break
            elif line.startswith("output_directory ="):
                subdir = f"Th{th}/Th{th}_U{u}"
                updated_directory = os.path.join(out_dir, subdir)
                lines[iline] = f"output_directory = {updated_directory}\n"

        fit_config_base = os.path.basename(fit_config)
        base_name, ext = os.path.splitext(fit_config_base)
        new_filename = f"{base_name}_Th{th}_U{u}{ext}"

        with open( str(configs_path + "/"+new_filename), "w") as file2:
            file2.writelines(lines)

    # Copy config files
    if event_config != "":
        os.system("cp " + str(event_config) + " " + str(configs_path) + "/" + os.path.basename(event_config) )
        event_config = configs_path + "/" + event_config.split("/")[-1]
    if pdf_config != "":
        os.system("cp " + str(pdf_config) + " " + str(configs_path) + "/" + os.path.basename(pdf_config) )
        pdf_config = configs_path + "/" + pdf_config.split("/")[-1]
    if syst_config != "":
        os.system("cp " + str(syst_config) + " " + str(configs_path) + "/" + os.path.basename(syst_config) )
        syst_config = configs_path + "/" + syst_config.split("/")[-1]
    if osc_config != "":
        os.system("cp " + str(osc_config) + " " + str(configs_path) + "/" + os.path.basename(osc_config) )
        osc_config = configs_path + "/" + osc_config.split("/")[-1]

    other_commands = f'sleep $[($RANDOM%{sleep_time+1})+1]s'

    # Write sh file
    out_macro_text = "#!/usr/bin/sh  \n" + \
                     "source " + env_file + "\n" + \
                     "source " + run_dir + "/env_my.sh " + "\n" + \
                     "cd " + str(run_dir) + "\n" + \
                     other_commands + "\n"
    for iFit in range(start_idx, end_idx):
        u = float(u_min) + iFit*(float(u_max)-float(u_min))/numvalsu
        u = "{:.3f}".format(u)

        fit_config_base = os.path.basename(fit_config)
        base_name, ext = os.path.splitext(fit_config_base)
        new_filename = f"{base_name}_Th{th}_U{u}{ext}"
        new_filename = configs_path + "/" + new_filename

        out_macro_text += f"{exec_path} {new_filename} {event_config} {pdf_config} {syst_config} {osc_config}\n"

    sh_filepath = "{0}/sh/".format(condor_path) + str(out_dir_base).replace("/", "") + '.sh'
    if not os.path.exists(os.path.dirname(sh_filepath)):
        os.makedirs(os.path.dirname(sh_filepath))
    with open(sh_filepath, "w") as sh_file:
        sh_file.write(out_macro_text)
    os.chmod(sh_filepath, 0o777)
    
    # Create submit file
    error_path = os.path.abspath('{0}/error'.format(condor_path))
    output_path = os.path.abspath('{0}/output'.format(condor_path))
    log_path = os.path.abspath('{0}/log'.format(condor_path))
    submit_path = os.path.abspath('{0}/submit'.format(condor_path))

    universe = "vanilla"
    notification = "never"
    n_rep = 1
    getenv = "True"

    submit_filepath = os.path.join(submit_path, out_dir_base)
    submit_filepath += ".submit"
    out_submit_text = \
        f"executable               = {sh_filepath}\n" \
        f"universe                 = {universe}\n" \
        f"output                   = {output_path}/{out_dir_base}.output\n" \
        f"error                    = {error_path}/{out_dir_base}.error\n" \
        f"log                      = {log_path}/{out_dir_base}.log\n" \
        f"notification             = {notification}\n" \
        f"priority                 = {priority}\n" \
        f"getenv                   = {getenv}\n" \
        f"allowed_execute_duration = {walltime} \n" \
        f"request_memory           = {mem} \n" \
        f"requirements             = (TARGET.CpuFamily == 25) \n" \
        f"queue {n_rep}\n"

    os.makedirs(os.path.dirname(submit_filepath), exist_ok=True)
    with open(submit_filepath, "w") as out_submit_file:
        out_submit_file.write(out_submit_text)

    command = f'condor_submit -batch-name \"{batch_name}\" {submit_filepath}'
    print ("executing job: " + command)
    os.system(command)

if __name__ == "__main__":

    parser = argparse.ArgumentParser("Submit fixed geo rate jobs to condor")
    parser.add_argument('exe', type=str, help='executable')
    parser.add_argument('out_dir', type=str, help='output directory')
    parser.add_argument('-r', "--run_dir", type=str, default=".", help='directory to run executable from')
    parser.add_argument('-e', "--env_file", type=str, default="", help='environment variable file to source')
    parser.add_argument('-f', "--fit_cfg", type=str, default="", help='fit config path')
    parser.add_argument('-i', "--event_cfg", type=str, default="", help='event config path')
    parser.add_argument('-p', "--pdf_cfg", type=str, default="", help='pdf config path')
    parser.add_argument('-s', "--syst_cfg", type=str, default="", help='syst config path')
    parser.add_argument('-o', "--osc_cfg", type=str, default="", help='osc grid config path')
    parser.add_argument('-d', "--numu", type=int, default=250, help='number of U points (individual fits) per job')
    parser.add_argument("-w", "--wall_time", type=int, default=86400, help="max runtime (seconds)")
    parser.add_argument("-m", "--mem", type=float, default=300, help="max memory (MB)")
    parser.add_argument("-j", "--job_name", type=str, default="", help='job name')
    args = parser.parse_args()

    exec_name = args.exe
    out_dir = check_dir(args.out_dir)
    base_name = out_dir.split("/")[-2]
    run_dir = args.run_dir
    env_file = args.env_file
    fit_config = run_dir + "/" + args.fit_cfg if args.fit_cfg != "" else ""
    event_config = run_dir + "/" + args.event_cfg if args.event_cfg != "" else ""
    pdf_config = run_dir + "/" + args.pdf_cfg if args.pdf_cfg != "" else ""
    syst_config = run_dir + "/" + args.syst_cfg if args.syst_cfg != "" else ""
    osc_config = run_dir + "/" + args.osc_cfg if args.osc_cfg != "" else ""
    u_chunk = args.numu
    walltime = args.wall_time
    mem = args.mem
    job_name = args.job_name

    u_min, u_max, th_min, th_max = read_fitcfg(fit_config)

    if job_name == "":
        job_name = base_name

    for iTh in range(numvalsth):
        th = float(th_min) + iTh*(float(th_max)-float(th_min))/numvalsth
        th = "{:.3f}".format(th)

        for start_idx in range(0, numvalsu, u_chunk):
            end_idx = min(start_idx + u_chunk, numvalsu)
            batch_name = f"{job_name}_Th{th}_U{start_idx}-{end_idx-1}"
            log_dir = check_dir("{0}/log/".format(out_dir))
            error_dir = check_dir("{0}/error/".format(out_dir))
            sh_dir = check_dir("{0}/sh/".format(out_dir))
            submit_dir = check_dir("{0}/submit/".format(out_dir))
            output_dir = check_dir("{0}/output/".format(out_dir))

            pycondor_submit(batch_name, exec_name, out_dir, run_dir, env_file,
                            fit_config, event_config, pdf_config, syst_config,
                            osc_config, walltime, mem, th,
                            start_idx, end_idx, sleep_time = 1, priority = 5)
