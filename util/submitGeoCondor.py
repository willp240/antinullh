import os
import argparse
from typing import List, Tuple, Optional

ratio_min = 0.1
ratio_max = 10.0
ratio_step = 0.1

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

def makeRatioList(rmin = 0.1, rmax = 10.0, rstep = 0.1):
    # Avoid float accumulation issues by stepping integers
    n = int(round((rmax - rmin) / rstep)) + 1
    for i in range(n):
        return [round(rmin + i * rstep, 10) for i in range(n)]

def ratioTag(r):
    # e.g. 0.1 -> "0p1", 10.0 -> "10p0"
    s = f"{r:.1f}"
    return s.replace(".", "p")

def parseNom(lines, block_name):
    """
    Find [block_name] and return float value of 'nom = ...' within that block.
    Stops at next [....] header.
    """
    in_block = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            in_block = (stripped == f"[{block_name}]")
            continue
        if in_block and stripped.startswith("nom"):
            parts = stripped.split("=", 1)
            if len(parts) == 2:
                return float(parts[1].strip())
    return None

def updateNom(lines, block_name, new_nom):
    """
    Update 'nom =' in a given [block_name] block. Returns modified lines.
    """
    in_block = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            in_block = (stripped == f"[{block_name}]")
            continue
        if in_block and stripped.startswith("nom"):
            lines[i] = f"nom = {new_nom}\n"
            break
    return lines

def updateDir(lines, out_dir, subdir):
    """
    Replace the top-level 'output_directory =' line in the fit cfg with out_dir/subdir.
    """
    updated_directory = os.path.join(out_dir, subdir)
    for i, line in enumerate(lines):
        if line.strip().startswith("output_directory"):
            lines[i] = f"output_directory = {updated_directory}\n"
            break
    return lines

def computeRates(total, ratio):
    # ratio r = U/Th, total S = U + Th
    # => Th = S/(1+r), U = r*S/(1+r)
    th = total / (1.0 + ratio)
    u = ratio * th
    return u, th

def pycondor_submit(job_name, exec_name, out_dir, run_dir, env_file,
                    fit_config, event_config, pdf_config, syst_config, osc_config,
                    walltime, mem, ratio, sleep_time = 1, priority= 5):
    """
    Submit one job that runs one fit with U/Th fixed to the given ratio (sum held constant).
    Writes a per-job fit cfg under out_dir/ratioX/cfg and output_directory under out_dir/ratioX.
    """
    print(job_name)
    batch_name = job_name.split('_', 1)[0]

    condor_path = f"{out_dir}/"
    exec_path = os.path.join(run_dir, "bin", exec_name)

    rtag = ratioTag(ratio)
    ratio_dir = f"ratio{rtag}"

    configs_path = os.path.abspath(f"{condor_path}/{ratio_dir}/cfg")
    check_dir(configs_path)

    with open(fit_config, "r") as f:
        base_lines = f.readlines()

    # Get totals from input cfg (PPO and bis-MSB blocks if present)
    u0  = parseNom(base_lines, "geonu_U")
    th0 = parseNom(base_lines, "geonu_Th")
    u20  = parseNom(base_lines, "geonu_U2")
    th20 = parseNom(base_lines, "geonu_Th2")

    if u0 is None or th0 is None:
        raise RuntimeError("Could not find nom values for [geonu_U] and/or [geonu_Th] in fit cfg.")

    total1 = u0 + th0
    new_u, new_th = computeRates(total1, ratio)

    # Work on a fresh copy of lines for this ratio
    lines = list(base_lines)

    # Update nom values for dataset 1
    lines = updateNom(lines, "geonu_U", new_u)
    lines = updateNom(lines, "geonu_Th", new_th)

    # Update dataset 2 if present
    if u20 is not None and th20 is not None:
        total2 = u20 + th20
        new_u2, new_th2 = computeRates(total2, ratio)
        lines = updateNom(lines, "geonu_U2", new_u2)
        lines = updateNom(lines, "geonu_Th2", new_th2)

    # Update output directory for this job
    lines = updateDir(lines, out_dir, ratio_dir)

    # Write new fit cfg filename
    fit_config_base = os.path.basename(fit_config)
    base_name, ext = os.path.splitext(fit_config_base)
    new_fit_filename = f"{base_name}_ratio{rtag}{ext}"
    new_fit_path = os.path.join(configs_path, new_fit_filename)

    with open(new_fit_path, "w") as f2:
        f2.writelines(lines)

    # Copy config files into cfg dir (same behaviour as your script)
    if event_config != "":
        os.system("cp " + str(event_config) + " " + str(configs_path) + "/" + os.path.basename(event_config))
        event_config = os.path.join(configs_path, os.path.basename(event_config))
    if pdf_config != "":
        os.system("cp " + str(pdf_config) + " " + str(configs_path) + "/" + os.path.basename(pdf_config))
        pdf_config = os.path.join(configs_path, os.path.basename(pdf_config))
    if syst_config != "":
        os.system("cp " + str(syst_config) + " " + str(configs_path) + "/" + os.path.basename(syst_config))
        syst_config = os.path.join(configs_path, os.path.basename(syst_config))
    if osc_config != "":
        os.system("cp " + str(osc_config) + " " + str(configs_path) + "/" + os.path.basename(osc_config))
        osc_config = os.path.join(configs_path, os.path.basename(osc_config))

    other_commands = f"sleep $[($RANDOM%{sleep_time+1})+1]s"

    # Write sh file (ONE fit)
    out_macro_text = (
        "#!/usr/bin/sh  \n"
        + "source " + env_file + "\n"
        + "source " + run_dir + "/env_my.sh " + "\n"
        + "cd " + str(run_dir) + "\n"
        + other_commands + "\n"
        + f"{exec_path} {new_fit_path} {event_config} {pdf_config} {syst_config} {osc_config}\n"
    )

    out_dir_base = f"{out_dir.split('/')[-2]}_{ratio_dir}"
    sh_filepath = f"{condor_path}/sh/" + str(out_dir_base).replace("/", "") + ".sh"
    os.makedirs(os.path.dirname(sh_filepath), exist_ok=True)
    with open(sh_filepath, "w") as sh_file:
        sh_file.write(out_macro_text)
    os.chmod(sh_filepath, 0o777)

    # Create submit file
    error_path = os.path.abspath(f"{condor_path}/error")
    output_path = os.path.abspath(f"{condor_path}/output")
    log_path = os.path.abspath(f"{condor_path}/log")
    submit_path = os.path.abspath(f"{condor_path}/submit")

    universe = "vanilla"
    notification = "never"
    getenv = "True"
    n_rep = 1

    submit_filepath = os.path.join(submit_path, out_dir_base + ".submit")
    out_submit_text = (
        f"executable               = {sh_filepath}\n"
        f"universe                 = {universe}\n"
        f"output                   = {output_path}/{out_dir_base}.output\n"
        f"error                    = {error_path}/{out_dir_base}.error\n"
        f"log                      = {log_path}/{out_dir_base}.log\n"
        f"notification             = {notification}\n"
        f"priority                 = {priority}\n"
        f"getenv                   = {getenv}\n"
        f"allowed_execute_duration = {walltime}\n"
        f"request_memory           = {mem}\n"
        f"queue {n_rep}\n"
    )

    os.makedirs(os.path.dirname(submit_filepath), exist_ok=True)
    with open(submit_filepath, "w") as out_submit_file:
        out_submit_file.write(out_submit_text)

    command = f'condor_submit -batch-name "{batch_name}" {submit_filepath}'
    print("executing job: " + command)
    os.system(command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Submit geo-U/Th ratio scan jobs to condor")
    parser.add_argument("exe", type=str, help="executable")
    parser.add_argument("out_dir", type=str, help="output directory")
    parser.add_argument("-r", "--run_dir", type=str, default=".", help="directory to run executable from")
    parser.add_argument("-e", "--env_file", type=str, default="", help="environment variable file to source")
    parser.add_argument("-f", "--fit_cfg", type=str, default="", help="fit config path")
    parser.add_argument("-i", "--event_cfg", type=str, default="", help="event config path")
    parser.add_argument("-p", "--pdf_cfg", type=str, default="", help="pdf config path")
    parser.add_argument("-s", "--syst_cfg", type=str, default="", help="syst config path")
    parser.add_argument("-o", "--osc_cfg", type=str, default="", help="osc grid config path (passed through to executable)")
    parser.add_argument("-n", "--numratios", type=int, default=1, help=" 'ratios per job' (default 1)")
    parser.add_argument("-w", "--wall_time", type=int, default=86400, help="max runtime (seconds)")
    parser.add_argument("-m", "--mem", type=float, default=300, help="max memory (MB)")
    parser.add_argument("-j", "--job_name", type=str, default="", help="job name")
    args = parser.parse_args()

    exec_name = args.exe
    out_dir = check_dir(args.out_dir)
    base_name = out_dir.split("/")[-2]
    run_dir = args.run_dir
    env_file = args.env_file
    fit_config = os.path.join(run_dir, args.fit_cfg) if args.fit_cfg != "" else ""
    event_config = os.path.join(run_dir, args.event_cfg) if args.event_cfg != "" else ""
    pdf_config = os.path.join(run_dir, args.pdf_cfg) if args.pdf_cfg != "" else ""
    syst_config = os.path.join(run_dir, args.syst_cfg) if args.syst_cfg != "" else ""
    osc_config = os.path.join(run_dir, args.osc_cfg) if args.osc_cfg != "" else ""
    walltime = args.wall_time
    mem = args.mem
    job_name = args.job_name
    if job_name == "":
        job_name = base_name

    ratios = makeRatioList(ratio_min, ratio_max, ratio_step)
    ratios_per_job = max(1, int(args.numratios))

    for r in ratios:
        rtag = ratioTag(r)
        batch_name = f"{job_name}_ratio{rtag}"
        log_dir = check_dir("{0}/log/".format(out_dir))
        error_dir = check_dir("{0}/error/".format(out_dir))
        sh_dir = check_dir("{0}/sh/".format(out_dir))
        submit_dir = check_dir("{0}/submit/".format(out_dir))
        output_dir = check_dir("{0}/output/".format(out_dir))
        pycondor_submit(batch_name, exec_name, out_dir, run_dir, env_file,
                        fit_config, event_config, pdf_config, syst_config, osc_config,
                        walltime, mem, r, sleep_time=1, priority=5)

