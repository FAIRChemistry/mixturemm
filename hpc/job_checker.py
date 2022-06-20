"""Script for checking on jobs.

MIT License

Copyright (c) 2021, Benjamin Schmitz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE."""

import os
import sys
import pandas as pd
import statistics

# get arguments from the command-line this file was started with
nvt_duration_ns = float(sys.argv[1])
nve_duration_ns = float(sys.argv[2])
if sys.argv[3] == 'True':
    verbose = True
else:
    verbose = False
energy_shift_tolerance_percent = float(sys.argv[4])
npt_duration_ns = float(sys.argv[5])
total_duration = nvt_duration_ns + nve_duration_ns
workdir = os.getcwd()
energy_warning_list = []

def try_convert_to_float(string):
    """Tries to convert a string into a float, if that is not possible returns None.

    Args:
        string (str): A string that should try to be converted to a float.

    Returns:
        float: The resulting float of the conversion.
    """
    try:
        return float(string)
    except:
        return None

def get_folder_strings_from_file(path_to_list_file):
    """Reads the folder paths on the high-performance computing center from a list text file with folder path strings.

    Args:
        path_to_list_file (str): Path to a list text file with folder path strings on the high-performance computing center.

    Returns:
        list: List with all replica folder paths.
    """

    with open(f'{path_to_list_file}') as file:
        return [line.rstrip() for line in file]

def get_run_duration(path_to_csv_file):
    """Extracts the run duration of a replica from the csv file obtained through the StateDataReporter of OpenMM.

    Args:
        path_to_csv_file (str): Path to the csv file obtained by the StateDataReporter of OpenMM.

    Returns:
        float: Run duration in nanoseconds.
    """

    replica_df = pd.read_csv(path_to_csv_file)
    first_line = replica_df.head(n=1)
    last_line = replica_df.tail(n=1)
    first_line.reset_index(inplace=True)
    last_line.reset_index(inplace=True)
    start_time_ = first_line.loc[0].at['#"Time (ps)"']
    progress_time_ = last_line.loc[0].at['#"Time (ps)"']
    start_time = float(start_time_)
    progress_time = float(progress_time_)
    run_duration = (progress_time - start_time) * 10**(-3)
    return run_duration

def get_avg_temperature(path_to_csv_file):
    """Extracts the average temperature of a replica from the csv file obtained through the StateDataReporter of OpenMM.

    Args:
        path_to_csv_file (str): Path to the csv file obtained by the StateDataReporter of OpenMM.

    Returns:
        float: Average temperature and standard deviation in Kelvin.
    """

    replica_df = pd.read_csv(path_to_csv_file)
    temperature_list__ = replica_df.loc[:, 'Temperature (K)'].to_list()
    temperature_list_ = [temp if isinstance(temp, (int, float)) else try_convert_to_float(temp) for temp in temperature_list__]
    temperature_list = [temp for temp in temperature_list_ if temp is not None]
    return statistics.mean(temperature_list), statistics.stdev(temperature_list)

def get_total_energy_shift(path_to_csv_file):
    """Extracts the total energy shift of a replica from the csv file obtained through the StateDataReporter of OpenMM.

    Args:
        path_to_csv_file (str): Path to the csv file obtained by the StateDataReporter of OpenMM.

    Returns:
        float: Total energy shift in percent.
    """

    replica_df = pd.read_csv(path_to_csv_file)
    total_energy_list__ = replica_df.loc[:, 'Total Energy (kJ/mole)'].to_list()
    total_energy_list_ = [ene if isinstance(ene, (int, float)) else try_convert_to_float(ene) for ene in total_energy_list__]
    total_energy_list = [ene for ene in total_energy_list_ if ene is not None]
    return abs((statistics.stdev(total_energy_list)/statistics.mean(total_energy_list))) * 100

def write_job_checks_header(path_to_job_checks):
    """Writes a header for the job_checks text file.

    Args:
        path_to_job_checks (str): Path to the job_checks text file.
    """

    with open(f'{path_to_job_checks}', 'w') as file:
        file.write('''
 _______     _        _______ _                 _          
(_______)   | |      (_______) |               | |         
     _  ___ | |__     _      | |__  _____  ____| |  _  ___ 
 _  | |/ _ \|  _ \   | |     |  _ \| ___ |/ ___) |_/ )/___)
| |_| | |_| | |_) )  | |_____| | | | ____( (___|  _ (|___ |
 \___/ \___/|____/    \______)_| |_|_____)\____)_| \_|___/ 
                                                           
===========================================================


''')

def write_job_checks_entry_system(folder_string, path_to_job_checks, npt_duration_ns):
    """Writes an entry for a system into the job_checks text file containing progress percentage.

    Args:
        folder_string (str): System folder path.
        path_to_job_checks (str): Path to the job_checks text file.
        npt_duration_ns (float): NpT equilibration duration in nanoseconds.
    """

    temperature = float(folder_string.split('/')[-1].split('-')[3])
    with open(f'{path_to_job_checks}', 'a') as file:
        if os.path.exists(f'{folder_string}/npt_equilibration/npt_equilibration_end.pdb'):
            file.write(f'{folder_string} density adjusted\n')
            file.write(f'System simulation at 100 %\n')
        else:
            file.write(f'{folder_string} not yet adjusted\n')
            if os.path.exists(f'{folder_string}/npt_equilibration/npt_equilibration.csv'):
                run_duration = get_run_duration(f'{folder_string}/npt_equilibration/npt_equilibration.csv')
                total_duration = npt_duration_ns + (temperature * 10**(-3))
                percentage_run = (run_duration / total_duration) * 100
                file.write(f'System simulation at {percentage_run:.2f} %\n')
            else:
                file.write(f'System simulation at 0 %\n')

def write_job_checks_entry_replica(folder_string, path_to_job_checks, nvt_duration_ns, total_duration, energy_shift_tolerance_percent, verbose):
    """Writes an entry for a replica into the job_checks text file containing progress percentage and optionally average temperature and total energy shift.

    Args:
        folder_string (str): Replica folder path.
        path_to_job_checks (str): Path to the job_checks text file.
        nvt_duration_ns (int/float): NVT equilibration duration in nanoseconds.
        total_duration (int/float): Added time of NVT equilibration and NVE production.
        energy_shift_tolerance_percent (float): Total energy shift tolerance in percent.
        verbose (bool): If true includes average temperature and total energy shift in the job_checker entry.
    """

    with open(f'{path_to_job_checks}', 'a') as file:
        if os.path.exists(f'{folder_string}/nve_production_end.pdb'):
            file.write(f'{folder_string} done\n')
            file.write(f'Replica simulation at 100 %\n')
            if verbose:
                avg_temp, stdev_temp = get_avg_temperature(f'{folder_string}/nve_production.csv')
                energy_shift = get_total_energy_shift(f'{folder_string}/nve_production.csv')
                if energy_shift > energy_shift_tolerance_percent:
                    energy_warning_list.append(f'{folder_string}: {energy_shift:.6f}')
                file.write(f'Replica temperature at {avg_temp:.2f} +/- {stdev_temp:.2f} K\n')
                file.write(f'Total energy shift at {energy_shift:.6f} %\n')
        else:
            file.write(f'{folder_string} not yet finished\n')
            if os.path.exists(f'{folder_string}/nve_production.csv'):
                run_duration_ = get_run_duration(f'{folder_string}/nve_production.csv')
                run_duration = run_duration_ + nvt_duration_ns
                percentage_run = (run_duration / total_duration) * 100
                file.write(f'Replica simulation at {percentage_run:.2f} %\n')
            elif os.path.exists(f'{folder_string}/nvt_equilibration.csv'):
                run_duration = get_run_duration(f'{folder_string}/nvt_equilibration.csv')
                percentage_run = (run_duration / total_duration) * 100
                file.write(f'Replica simulation at {percentage_run:.2f} %\n')
            else:
                file.write(f'Replica simulation at 0 %\n')

def write_job_checks_footer(path_to_job_checks):
    """Writes a concluding footer into the job_checks text file containing the number of planned, done and not yet finished replicas and energy shift warnings if the total energy shift tolerance has been exceeded.

    Args:
        path_to_job_checks (str): Path to the job_checks text file.
    """

    with open(f'{path_to_job_checks}', 'r') as file:
        data = file.read()
    count_done_systems = data.count('density adjusted')
    count_not_done_systems = data.count('not yet adjusted')
    total_systems = count_done_systems + count_not_done_systems
    count_done_replicas = data.count('done')
    count_not_done_replicas = data.count('not yet finished')
    total_replicas = count_done_replicas + count_not_done_replicas
    with open (f'{path_to_job_checks}', 'a') as file:
        file.write('''
=============================================
  ______                                     
 / _____)                                    
( (____  _   _ ____  ____  _____  ____ _   _ 
 \____ \| | | |    \|    \(____ |/ ___) | | |
 _____) ) |_| | | | | | | / ___ | |   | |_| |
(______/|____/|_|_|_|_|_|_\_____|_|    \__  |
                                      (____/ 
                                      
''')
        file.write('============================================================================================================================================\n')
        file.write(f'Total systems: {total_systems}\n')
        file.write(f'Density adjusted systems: {count_done_systems}\n')
        file.write(f'Not yet adjusted systems: {count_not_done_systems}\n')
        file.write(f'Total replicas: {total_replicas}\n')
        file.write(f'Done replicas: {count_done_replicas}\n')
        file.write(f'Not yet finished replicas: {count_not_done_replicas}\n')
        if not energy_warning_list:
            file.write(f'No energy shift warnings.\n')
        else:
            file.write(f'Energy shift warnings:\n')
            for warning_string in energy_warning_list:
                file.write(f'{warning_string}\n')
        file.write('============================================================================================================================================\n')

# This is the part that runs
write_job_checks_header(f'{workdir}/job_checks.txt')
system_folder_strings = get_folder_strings_from_file(f'{workdir}/system_list.txt')
replica_folder_strings = get_folder_strings_from_file(f'{workdir}/replica_list.txt')
for folder_string in system_folder_strings:
    write_job_checks_entry_system(folder_string, f'{workdir}/job_checks.txt', npt_duration_ns)
for folder_string in replica_folder_strings:
    if folder_string == 'finalize':
        write_job_checks_footer(f'{workdir}/job_checks.txt')
    else:
        write_job_checks_entry_replica(folder_string, f'{workdir}/job_checks.txt', nvt_duration_ns, total_duration, energy_shift_tolerance_percent, verbose)