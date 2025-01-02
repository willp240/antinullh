#!/usr/bin/python

#############################################
# Calculate the livetime for an input run list.
# The following information is output:
# Total livetime
# Average livetime
# Livetime info per run
############################################
# Information needed to run this script:
#   A run list, usually the golden run list is used.
#   A rat environment.
###############################################

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals

from datetime import datetime, timedelta
try:
    input = raw_input  # Python 2
    import pytz as timezone
except NameError:
    from datetime import timezone
    pass  # Python 3
import argparse
import os
import sys
import glob
import numpy as np
from rat import ratdb
import rat

# contact ratdb
RATDB_USER='snoplus'
RATDB_PASS="dontestopmenow"
RATDB_SERVER='pgsql.snopl.us'
RATDB_PORT=5400

livetime_dic = {
'muon': 0.0,
'raw_livetime': 0.0,
'missedmuon': 0.0,
'livetime': 0.0,
'overlap': 0.0,
}

data_files = []

#Format KEY=run, VAL=2D array of vetos for run with data
#Veto data array in the form: GTID, Nhit, clockCount50, clockCount10, UTDays, UTSecs, UTNSecs, Veto length (in ns)
veto_dic = {}

do_file_checks = True
do_muon_veto = True
do_bespoke_veto = True

def make_veto_dic(veto_dir,run_number):
    #Bespoke veto outputs are assumed to be contained in text files with
    #ending pattern below. Avoid putting duplicate veto processing
    #for a single run in the same directory, or any other
    #files besides veto outputs for that matter.
    veto_files = glob.glob(veto_dir+str(run_number)+'*.txt')

    for vfi in veto_files:
        veto_data = np.loadtxt(vfi,dtype=float,delimiter=' ', skiprows=1)
        print("veto_datandim ",veto_data.ndim)
        if veto_data.ndim == 1 and len(veto_data) != 0:
            veto_dic[veto_data[0]] = [veto_data[1:]]
        else:
            for veto in veto_data:
                if veto[0] not in veto_dic:
                    veto_dic[veto[0]] = [veto[1:]]
                else:
                    veto_dic[veto[0]].append(veto[1:])

def is_table_in_ratdb(run_number):
    '''
    Check if the run table exists for a specified run.
    If it does, calculate the various livetime values and write to a file.
    '''
    try:
        ratdbLink = ratdb.RATDBConnector('postgres://'+RATDB_USER+':'+RATDB_PASS+'@'+RATDB_SERVER+'/ratdb')
        table = ratdbLink.fetch(obj_type='RUN', run=int(run_number))
        run_table = table[0]['data']

        subfile_list = run_table['subfile_list']        
        start_day = run_table['start_day']
        start_sec = run_table['start_sec']
        start_nsc = run_table['start_nsc']
        stop_day = run_table['stop_day']
        stop_sec = run_table['stop_sec']
        stop_nsc = run_table['stop_nsc']

        missing_subruns = []
        if do_file_checks:
            run_str = 'r{:010d}'.format(run_number)
            subrun_list = [run_str + '_s' + subfi.split('_')[-1].split('.zdab')[0] for subfi in subfile_list]# i.e. 300000_s_000
            for subrun in subrun_list:
                if subrun not in data_files:
                    print('Run ' + str(run_number) + ': Subrun not found: ' + subrun)
                    missing_subruns.append(subrun)

        #SNO+ time counts from 1 Jan 2010 at 0:00 GMT (UTC 0)
        epoch_time = datetime(2010, 1, 1, tzinfo=timezone.utc)
        start_time = epoch_time + timedelta(days=start_day, seconds=start_sec, microseconds=start_nsc*1e-3)
        stop_time = epoch_time + timedelta(days=stop_day, seconds=stop_sec, microseconds=stop_nsc*1e-3)
        raw_run_len = (stop_time - start_time).total_seconds()/60/60/24
        
        time_correction = timedelta()
        
        #Do veto corrections (muon + high nhit)
        veto_times = []

        #Gathers muon veto times from LAST_MUON and TPMUONFOLLOWER tables
        if do_muon_veto:
            #Currently only consider short time, as long is unused as of Jan 2024
            #Will need to add new veto logic for long if revived
            muon_short_length = timedelta(seconds=20)

            table_lastmuon = ratdbLink.fetch(obj_type='LAST_MUON', run=int(run_number))
            try:
                lastmuon_table = table_lastmuon[0]['data']
                if lastmuon_table['day_last_muon'] != -1:
                    lm_day = lastmuon_table['day_last_muon']
                    lm_sec = lastmuon_table['sec_last_muon']
                    lm_ns = lastmuon_table['ns_last_muon']
                    lastmuon_start_time = epoch_time + timedelta(days=lm_day, seconds=lm_sec, microseconds=lm_ns*1e-3)
                else:
                    lastmuon_start_time = epoch_time
                lastmuon_end_time = lastmuon_start_time + muon_short_length
                if lastmuon_end_time > start_time:
                    veto_times.append((lastmuon_start_time, lastmuon_end_time))
            except IndexError:
                lastmuon_start_time = epoch_time

            table_tpmuon = ratdbLink.fetch(obj_type='TPMUONFOLLOWER', run=int(run_number))
            tpmuon_table = table_tpmuon[0]['data']
            tpmuon_days = tpmuon_table['day_muons']
            tpmuon_secs = tpmuon_table['sec_muons']
            tpmuon_ns = tpmuon_table['ns_muons']

            muon_veto_start_time_prev = lastmuon_start_time
            if tpmuon_days != -1:
                if len(tpmuon_days) != len(tpmuon_secs) or len(tpmuon_days) != len(tpmuon_ns):
                    print('Run ' + str(run_number) + ': Muon time arrays not the same length. Weird!!!')
                for day, sec, ns in zip(tpmuon_days, tpmuon_secs, tpmuon_ns):
                    muon_veto_start_time = epoch_time + timedelta(days=day, seconds=sec, microseconds=ns*1e-3)
                    muon_veto_end_time = muon_veto_start_time + muon_short_length
                    #If the muon times are not in order in the table, something weird has happened
                    if muon_veto_start_time < muon_veto_start_time_prev:
                        print('Run ' + str(run_number) + ': Muon veto start time decreased compared to previous. Weird!!!')
                        print(muon_veto_start_time_prev, muon_veto_start_time)
                    veto_times.append((muon_veto_start_time, muon_veto_end_time))
                    muon_veto_start_time_prev = muon_veto_start_time

        #Gathers bespoke veto times from specified output files
        if do_bespoke_veto:
            veto_start_time_prev = epoch_time
            table_veto_prev = veto_dic[run_number-1] if run_number-1 in veto_dic else []
            if table_veto_prev == []:
                print('Run ' + str(run_number) + ': Previous veto table empty!')

            for veto in table_veto_prev:
                veto_start_time = epoch_time + timedelta(days=veto[4], seconds=veto[5], microseconds=veto[6]*1e-3)
                veto_end_time = veto_start_time + timedelta(microseconds=veto[7]*1e-3)
                if veto[4] == 0 and veto[5] == 0 and veto[6] == 0:
                    print('Run ' + str(run_number) + ': Veto time not valid! Continuing!')
                    continue
                if veto_start_time > start_time:
                    print('Run ' + str(run_number) + ': Previous run veto time after start time of current run. Weird!!!')
                if veto_start_time > stop_time:
                    print('Run ' + str(run_number) + ': Previous run veto start time is after stop time of current run. Weird!!!')
                if veto_start_time < veto_start_time_prev:
                    print('Run ' + str(run_number) + ': Previous run veto start time decreased compared to previous. Weird!!!')
                if veto_end_time > start_time:
                    veto_times.append((veto_start_time, veto_end_time))
                    veto_start_time_prev = veto_start_time

            
            table_veto = veto_dic[run_number] if run_number in veto_dic else []
            if table_veto == []:
                print('Run ' + str(run_number) + ': Veto table empty!')

            for veto in table_veto:
                veto_start_time = epoch_time + timedelta(days=veto[4], seconds=veto[5], microseconds=veto[6]*1e-3)
                veto_end_time = veto_start_time + timedelta(microseconds=veto[7]*1e-3)
                if veto[4] == 0 and veto[5] == 0 and veto[6] == 0:
                    print('Run ' + str(run_number) + ': Veto time not valid! Continuing!')
                    continue
                if veto_start_time > stop_time:
                    print('Run ' + str(run_number) + ': Veto start time is after stop time of current run. Weird!!!')
                if veto_start_time < veto_start_time_prev:
                    print('Run ' + str(run_number) + ': Veto start time decreased compared to previous. Weird!!!')
                veto_times.append((veto_start_time, veto_end_time))
                veto_start_time_prev = veto_start_time

        #Sort the vetos so that they are in time order
        #Assumes that all vetos should be treated equally
        #Can add hiearchical logic later if needed
        veto_times = sorted(veto_times)

        if len(veto_times) != 0:
            veto_chain_start_time = veto_times[0][0] if veto_times[0][0] > start_time else start_time
            veto_chain_end_time = veto_times[0][1] if veto_times[0][1] < stop_time else stop_time
            for veto in veto_times[1:]:
                #The proposed next veto in the chain starts at the next veto time (vetoes are sorted in ascending order)
                #and ends either at the next veto end or the end of the run, whichever is sooner
                next_veto_start = veto[0]
                next_veto_end = veto[1] if veto[1] < stop_time else stop_time

                #First, check whether the proposed start time is less than the end of the current chain
                #and whether the proposed end time is later than the end of the current chain
                if next_veto_start < veto_chain_end_time and next_veto_end >= veto_chain_end_time:
                    #If so, the new chain end time is the next veto end (start remains the same)
                    veto_chain_end_time = next_veto_end
                #Second, check whether the proposed start time is less than the end of the current chain
                #and whether the proposed end time is sooner than the end of the current chain
                elif next_veto_start < veto_chain_end_time and next_veto_end < veto_chain_end_time:
                    #If so, this veto is entirely within the chain already, so we just continue
                    continue
                #Last, this means the next veto start is after the current chain ends
                #So we add the time correction and start a new chain starting with this veto
                else:
                    time_correction += (veto_chain_end_time - veto_chain_start_time)
                    veto_chain_start_time = next_veto_start
                    veto_chain_end_time = next_veto_end
                #If the veto chain is ending equal to or past the run end time,
                #then no more deadtime can be accrued and we end the loop
                if veto_chain_end_time >= stop_time:
                    veto_chain_end_time = stop_time
                    break

            #Need to close out loop since last chain will not be added
            time_correction += (veto_chain_end_time - veto_chain_start_time)


        time_correction_len = time_correction.total_seconds()/60/60/24
        corrected_run_len = raw_run_len - time_correction_len

        livetime_dic['raw_livetime'] += raw_run_len
        livetime_dic['livetime'] += corrected_run_len
        writestring = str(run_number) + ' ' + str(raw_run_len) + ' ' + str(corrected_run_len)
        if raw_run_len < corrected_run_len:
            print('Run ' + str(run_number) + ': Raw livetime less than corrected livetime. Weird!!!')
        return True, writestring, missing_subruns
    except ValueError:
        return False, 'Table Not Found!', []

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i',type = str,
                        help='Input filelist containing line separated run numbers.')
    parser.add_argument('--data_dir', '-d', type = str, default = os.getcwd(),
                        help='Directory containing data files (assumed ntuple, default is current directory).')
    parser.add_argument('--lower_range', '-l', type = int, default = 0,
                        help='Lowest run of range to calculate livetime.')
    parser.add_argument('--upper_range', '-u', type = int, default = 9999999,
                        help='Uppermost run of range to calculate livetime.')
    parser.add_argument('--tag', '-t', type = str, default = 'defaultlivetime',
                        help='Name of the livetime calculation trial defined by user.')
    parser.add_argument('--output_dir', '-o', type = str, default = os.getcwd(),
                        help='Directory to place output files in, by default the current directory.')
    parser.add_argument('--veto_dir', type = str, default = os.getcwd(),
                        help='Directory to search for bespoke veto text files, by default the current directory.')
    parser.add_argument('--turn_off_file_checks', action='store_true',
                        help='Whether to turn off checks for missing files or subruns and multiple passes of same file (default is ON)')
    parser.add_argument('--turn_off_muon_veto', action='store_true',
                        help='Whether to turn off deadtime calculation from muon veto (default is ON)')
    parser.add_argument('--turn_off_bespoke_veto', action='store_true',
                        help='Whether to turn off deadtime calculation from bespoke veto (default is ON)')
    args = parser.parse_args()

    do_file_checks = not args.turn_off_file_checks
    if do_file_checks:
        print('Performing file checks.')
    else:
        print('Not performing file checks.')
    do_muon_veto = not args.turn_off_muon_veto
    if do_muon_veto:
        print('Including muon veto in livetime calculation.')
    else:
        print('Not including muon veto in livetime calculation.')
    do_bespoke_veto = not args.turn_off_bespoke_veto
    if do_bespoke_veto:
        print('Including bespoke veto in livetime calculation.')
    else:
        print('Not including bespoke veto in livetime calculation.')

    if not RATDB_PASS:
        print('No RATDB password specified. Please fill in.')
        sys.exit()

    # Create a new directory to write out the livetime information
    outdir = args.output_dir + '/' + args.tag + '/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    #Grab data files in appropriate directory and truncate off extra info
    data_files = glob.glob(args.data_dir + 'Analysis*.ntuple.root')
    data_files = ['r0'+fi.split('.ntuple.root')[0].split('_r0')[-1] for fi in data_files]

   

    # Get the run list from the input file
    with open(args.input, 'r') as fruns:
        run_list = [int(line) for line in fruns]

    # Subset of the run list using the option run-range
    lower_run = args.lower_range
    upper_run = args.upper_range
    run_array = np.array(run_list)
    run_array = run_array[np.logical_and(run_array>=lower_run, run_array<=upper_run)]
    total_run = run_array.flatten().shape[0]

    # Open livetime information files
    fmisrun = open(outdir + str(run_array[0])+'missing_runs_' + str(args.tag) + '.txt','w')
    fmissub = open(outdir + str(run_array[0])+'missing_subruns_' + str(args.tag) + '.txt','w')
    fmulpass = open(outdir +str(run_array[0])+ 'multiple_passes_' + str(args.tag) + '.txt','w')
    ftot = open(outdir +str(run_array[0])+ 'livetime_total_'   + str(args.tag) + '.txt','w')
    favg = open(outdir + str(run_array[0])+'livetime_average_' + str(args.tag) + '.txt','w')
    frun = open(outdir + str(run_array[0])+'livetime_per_run_' + str(args.tag) + '.txt','w')

    if do_file_checks:
        #Check that there are files for every run in the list and that there are no overlapping passes
        for i, run in enumerate(run_array):
            run_files = [fi for fi in data_files if 'r0000'+str(run) in fi]
            #Check whether files exist for run
            if len(run_files) == 0:
                print('No files for run: ' + str(run))
                fmisrun.write(str(run) + '\n')

            #Check no overlapping passes
            subruns = set([int(fi.split('_s')[-1].split('_p')[0]) for fi in run_files])
            pass_files = [[fi for fi in run_files if 'r0000'+str(run)+'_s'+'{:03d}'.format(subrun) in fi] for subrun in subruns]
            multi_pass = []
            pass_file_lens = [(print('Multiple passes: ' + str(files)), multi_pass.append(files)) for files in pass_files if len(files) > 1]
            [fmulpass.write(' '.join(files) + '\n') for files in multi_pass]
    data_files = [fi.split('_p')[0] for fi in data_files]

    #Make veto dictionary from offline files
    veto_dir = args.veto_dir if args.veto_dir[-1] == '/' else args.veto_dir + '/'
    make_veto_dic(args.veto_dir,run_array[0])

    # Loop over the runs and get the information from the livetime table
    print('Calculating livetime for run-list...')
    print('Run Livetime Fractional Livetime')
    print('--------------------------------')
    for i, run in enumerate(run_array):
        table_flag, writestring, missing_subruns = is_table_in_ratdb(run)
        print(i, run, writestring, missing_subruns)

        if len(missing_subruns) != 0:
            fmissub.write('\n'.join(missing_subruns) + '\n')
        if table_flag:
            frun.write(writestring + '\n')
        else:
            fmisrun.write(str(run) + '\n')

    for key, val in livetime_dic.items():
        ftot.write(str(key) + ': ' + str(val) + '\n')
        favg.write(str(key) + ': ' + str(val / total_run) + '\n')

    fmissub.close()
    fmulpass.close()
    fmisrun.close()
    frun.close()
    ftot.close()
    favg.close()
    print('Completed livetime calculations')