# aims to change 
# 1. fit_config: fake_data_val in either energy_scale or energy_conv
# 2. copying rest of the configs
# 3. execute command line  i.e. python util/submitFixedOscJobs.py fixedosc_fit 
#                           /data/snoplus/weiiiiiii/antinuFit/run_selection/2p2energyshift_3percent 
#                           -r /home/huangp/antinullh/ -e /data/snoplus/software/snocave_Alma9/set_env.sh 
#                           -f cfg/2p2energyshift_3percent/fit_config.ini -i cfg/2p2energyshift_3percent/event_config.ini 
#                           -p cfg/2p2energyshift_3percent/pdf_config.ini -s cfg/2p2energyshift_3percent/syst_config.ini 
#                           -o cfg/2p2energyshift_3percent/oscgrid_config.ini
# Example command: python3 util/makeconfig_forrunselection.py  
#                  -s /home/huangp/antinullh/cfg/2p2energyconv_5plus 
#                  -newdir /home/huangp/antinullh/cfg/2p2energyscale_1p5plus 
#                  -o /data/snoplus/weiiiiiii/antinuFit/run_selection/2p2energyscale_1p5plus 
#                  -scale 1.015 -conv 0.045


import os
import shutil
import configparser
import subprocess
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-s", help="source dir to copy configs", type=str)
parser.add_argument("-newdir", help="path to put new config", type=str)
parser.add_argument("-o", help="output directory of fits", type=str)
parser.add_argument("-scale", help="fake_data_val ofenergy scale ", type=float)
parser.add_argument("-conv", help="fake_data_val of energy convolution", type=float)
args = parser.parse_args()
source_directory = f"{args.s}"
destination_directory = args.newdir
basedirname = os.path.basename(args.newdir)
if os.path.exists(destination_directory):
    shutil.rmtree(destination_directory)
shutil.copytree(source_directory, destination_directory)
print(f"Copied directory {source_directory} to {destination_directory}")


# Initialize parser
config = configparser.ConfigParser()

# Read the config file
config.read(f"{args.newdir}/fit_config.ini")  # Replace with the actual file path


# Update values
config.set("summary", "output_directory", str(args.o))
config.set("energy_conv", "fake_data_val", str(args.conv))
config.set("energy_scale", "fake_data_val", str(args.scale))

# Save changes back to the file
with open(f"{args.newdir}/fit_config.ini", "w") as configfile:
    config.write(configfile)

print("Updated fake_data_val in [energy_conv] and [energy_scale].")


command = f"python util/submitFixedOscJobs.py fixedosc_fit {args.o} " \
          "-r /home/huangp/antinullh/ " \
          "-e /data/snoplus/software/snocave_Alma9/set_env.sh " \
          f"-f cfg/{basedirname}/fit_config.ini " \
          f"-i cfg/{basedirname}/event_config.ini " \
          f"-p cfg/{basedirname}/pdf_config.ini " \
          f"-s cfg/{basedirname}/syst_config.ini " \
          f"-o cfg/{basedirname}/oscgrid_config.ini"
print(command)
# Run the command
os.system(command)
print("Command is now running")
