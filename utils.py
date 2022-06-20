"""Module with useful functions for preparing a project.

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
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def directory_maker(path) -> str:
   """A function that creates directories.

   Args:
      path (str): Path to the directory that should be created.
   """

   if not os.path.exists(path):
      os.makedirs(path)
   return path
   
def add_chainsubmitter_script(path) -> None:
   """Adds a chainsubmitter.sh file to the selected folder that can be used to submit chain jobs on a HPC using SLURM.
   """

   with open(f'{path}/chainsubmitter.sh', 'w') as file:
      file.write(
'''#!/bin/bash
######################################################
##        submitter script for chain jobs           ##
######################################################

## Define maximum number of jobs via positional parameter 1, default is 5
max_nojob=${1:-5}

## Define location of the job scripts 
## as list of strings?
chain_link_job=${2:-${PWD}/shscript.sh}

## Define type of dependency via positional parameter 2, default is 'afterok'
dep_type="${3:-afternotok}"

## Define the queue you want to use
queue=${4:-gpu_4}

myloop_counter=1

while [ ${myloop_counter} -le ${max_nojob} ] ; do

## Differ msub_opt depending on chain link number
   if [ ${myloop_counter} -eq 1 ] ; then
      slurm_opt=""
   else
      slurm_opt="-d ${dep_type}:${jobID}"
   fi

## Print current iteration number and sbatch command
   echo "Chain job iteration = ${myloop_counter}"
   echo "sbatch --export=myloop_counter=${myloop_counter} ${slurm_opt} ${chain_link_job}"
   ## Store job ID for next iteration by storing output of sbatch command with empty lines
   jobID=$(sbatch -p ${queue} --export=ALL,myloop_counter=${myloop_counter} ${slurm_opt} ${chain_link_job} 2>&1 | sed 's/[S,a-z]* //g')

## Check if ERROR occured
   if [[ "${jobID}" =~ "ERROR" ]] ; then
      echo "   -> submission failed!" ; exit 1
   else
      echo "   -> job number = ${jobID}"
   fi

## Increase counter
   let myloop_counter+=1

done'''
)

def hpc_submission_header(file_name, 
                          process_name, 
                          hpc_workspace, 
                          environment_name,
                          scheduler,
                          partition,
                          number_of_threads,
                          number_of_gpus,
                          max_runtime_hh_mm_ss,
                          conda_module) -> None:
   """Writes the header of a hpc submission script.

   Args:
       file_name (str): The name of the hpc submission script.
       process_name (str): The process name that the queued job should get. 
       hpc_workspace (str): The path to the hpc workspace.
       environment_name (str): The name of the anaconda environment containing the dependencies.
       scheduler (str): The queue scheduler.
       partition (str): The partition the job should be submitted to.
       number_of_threads (int): The number of requested threads.
       number_of_gpus (int): The number of GPUs requested on a node.
       max_runtime_hh_mm_ss (str): Maximum runtime for a job.
       conda_module (bool): Adapts the header to whether there is a conda module available on the hpc or not.
   """
    
   with open(f'{file_name}', 'w', newline='\n') as file:
      file.write(
f'''#!/bin/bash
#{scheduler} -J {process_name}
''')
      if partition == 'accelerated':
         file.write(
f'''#{scheduler} --gres=gpu:{number_of_gpus}
''')
      file.write(
f'''#{scheduler} --ntasks={number_of_threads}
#{scheduler} --time={max_runtime_hh_mm_ss}
#{scheduler} --partition={partition}

''')
      if not conda_module:
         file.write(
f'''set +eu

module purge
module load compiler/pgi/2020
module load devel/cuda/11.0
source {hpc_workspace}/conda/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"
conda activate {environment_name}

set -eu

''')
      if conda_module:
         file.write(
f'''set +eu

module purge
module load compiler/pgi/2020
module load devel/cuda/11.0
module load devel/miniconda
eval "$(conda shell.bash hook)"
conda activate {environment_name}

set -eu

''')
      file.write(
f'''echo "-----------------------------------------------------------------------"
echo {process_name}
echo $(date -u) "Job was started"
echo "-----------------------------------------------------------------------"

''')

def merge_project_descriptions(path_to_old_description, 
                               new_json_structure) -> None:
   """Merges two project descriptions

   Args:
       path_to_old_description (str): The path to the already existent project description.
       new_json_structure (dict): The newly created json structure that should be merged with the existent one.
   """

   with open(f'{path_to_old_description}', 'r+', encoding='utf-8') as f:
      json_base = json.load(f)
      for key, value in json_base['project'].items():
         if key in new_json_structure['project']:
            for inner_key, inner_value in value.items():
                  if inner_key in new_json_structure['project'][f'{key}']:
                     if type(inner_value) is not list:
                        if inner_value != new_json_structure['project'][f'{key}'][f'{inner_key}']:
                              old = inner_value
                              new = new_json_structure['project'][f'{key}'][f'{inner_key}']
                              json_base['project'][f'{key}'][f'{inner_key}'] = [old, new]
                     else:
                        if type(new_json_structure['project'][f'{key}'][f'{inner_key}']) is not list:
                              if new_json_structure['project'][f'{key}'][f'{inner_key}'] not in inner_value:
                                 json_base['project'][f'{key}'][f'{inner_key}'].append(new_json_structure['project'][f'{key}'][f'{inner_key}'])
                        else:
                              inner_value.extend(x for x in new_json_structure['project'][f'{key}'][f'{inner_key}'] if x not in inner_value)
                              json_base['project'][f'{key}'][f'{inner_key}'] = inner_value
      f.seek(0)
      json.dump(json_base, f, ensure_ascii=False, indent=4)

def jsonize_project(Project) -> dict:
   """Creates the JSON structure of the project by putting project attributes in a dictionary.

   Args:
       Project (obj): The project holds all parameters and functions to do a complete systematic simulation of different systems in openmm.

   Returns:
       dict: The JSON structure of the project description.
   """

   json_structure = {'workflow_info': {'general': {'package_name': 'mixturemm',
                                                   'simulation_engine': 'OpenMM',
                                                   'analysis_package': 'MDAnalysis'},
                                       'simulation': {'ensemble_order': '---(slow heating)---> |NpT| ---(rescale boxes and assign velocities)---> |NVT| ------> |NVE|',
                                                      'pressure_and_temperature_control': {'NpT': 'Langevin thermostat, Monte Carlo barostat',
                                                                                           'NVT': 'Langevin thermostat',
                                                                                           'NVE': 'no pressure and temperature control'},
                                                      'constraints': 'hydrogen bonds',
                                                      'nonbonded forces algorithm': 'particle mesh ewald'}},
                     'project': {'simulation_parameters': {'total_number_of_molecules': Project.total_number_molecules,
                                                           'initial_box_side_lengths': Project.init_box_side_length,
                                                           'water_mole_fractions': Project.chi_water_s,
                                                           'temperature_s_in_kelvin': Project.temperature_s,
                                                           'pressure_s_in_bar': Project.npt_equilibration_pressure_s,
                                                           'number_of_replicas': Project.replica_count},
                                 'hardware_settings': {'platform': Project.simulation_platform,
                	                                     'properties': Project.simulation_properties},
                                 'simulation_settings': {'pme_error_tolerance': Project.pme_error_tolerance,
                                                         'constraint_tolerance': Project.constraint_tolerance,
                                                         'cutoff_distance_in_nm': Project.cutoff_distance_nm,
                                                         'switching_function_starting_distance_in_nm': Project.cutoff_switch_distance_nm},
                                 'npt_equilibration': {'pressure_coupling_frequency': Project.npt_equilibration_pressure_coupling_frequency,
                                                       'temperature_coupling_frequency': Project.npt_equilibration_temperature_coupling_frequency,
                                                       'timestep_in_fs': Project.npt_equilibration_timestep_fs,
                                                       'duration_in_ns': Project.npt_equilibration_duration_ns,
                                                       'state_data_reporting_frequency': Project.reporting_frequency_state_npt_equilibration},
                                 'nvt_equilibration': {'temperature_coupling_frequency': Project.nvt_equilibration_temperature_coupling_frequency,
                	                                     'timestep_in_fs': Project.nvt_equilibration_timestep_fs,
                                                       'duration_in_ns': Project.nvt_equilibration_duration_ns,
                                                       'state_data_reporting_frequency': Project.reporting_frequency_state_nvt_equilibration},
                                 'nve_production': {'timestep_in_fs': Project.nve_production_timestep_fs,
                                                    'duration_in_ns': Project.nve_production_duration_ns,
                                                    'state_data_reporting_frequency': Project.reporting_frequency_state_nve_production,
                                                    'unwrapped_trajectory_reporting_frequency': Project.reporting_frequency_coordinates_unwrapped,
                                                    'wrapped_trajectory_reporting_frequency': Project.reporting_frequency_coordinates_wrapped}}}
   return json_structure

def jsonize_molecule(Molecule) -> dict:
   """Creates the JSON structure of the molecule by putting project attributes in a dictionary.

   Args:
       Molecule (obj): Contains information on the molecule that is used for packing with packmol and to identify it.

   Returns:
       dict: The JSON structure of the project description.
   """

   json_structure = {'number_of_atoms': Molecule.number_of_atoms,
                     'abbreviation': Molecule.abbreviation,
                     'smiles_code': Molecule.smiles,
                     'inchi_key': Molecule.inchi,
                     'molar_mass': Molecule.molar_mass,
                     'used_as_water': Molecule.use_as_water}
   return json_structure

def jsonize_forcefield(Forcefield) -> dict:
   """Creates the JSON structure of the force field by putting project attributes in a dictionary.

   Args:
       Forcefield (obj): Contains information on the force field.

   Returns:
       dict: The JSON structure of the project description.
   """

   json_structure = {'built_in': Forcefield.built_in,
                     'added_elements': Forcefield.elements}
   return json_structure

def add_json_entry(path_to_json, 
                   entry, 
                   key=None, 
                   subkey=None, 
                   subsubkey=None,
                   subsubsubkey=None) -> None:
   """Adds an entry to a JSON file at the specified level defined by the given keys.

   Args:
       path_to_json (str): The path to the JSON file to that the entry should be added.
       entry (str/int/list/dict): The entry that should be added.
       key (str, optional): The first level key of the JSON. Defaults to None.
       subkey (str, optional): The second level key of the JSON. Defaults to None.
       subsubkey (str, optional): The third level key of the JSON. Defaults to None.
       subsubsubkey (str, optional): The fourth level key of the JSON. Defaults to None.
   """

   condition = sum(arg is not None for arg in [key, subkey, subsubkey, subsubsubkey])
   with open(f'{path_to_json}', 'r+', encoding='utf-8') as f:
      j_file = json.load(f)
      if condition == 1:
         j_file[key] = entry
      elif condition == 2:
         if key not in j_file:
            j_file[key] = {}
         j_file[key][subkey] = entry
      elif condition == 3:
         if key not in j_file:
            j_file[key] = {}
         if subkey not in j_file[key]:
            j_file[key][subkey] = {}
         j_file[key][subkey][subsubkey] = entry
      elif condition == 4:
         if key not in j_file:
            j_file[key] = {}
         if subkey not in j_file[key]:
            j_file[key][subkey] = {}
         if subsubkey not in j_file[key][subkey]:
            j_file[key][subkey][subsubkey] = {}
         j_file[key][subkey][subsubkey][subsubsubkey] = entry
      f.seek(0)
      json.dump(j_file, f, ensure_ascii=False, indent=4)

def bash_arg_prepper(*args) -> tuple:
   """Prepares the command line arguments for the HPC bash files.

   Args:
       *args: Variable length argument list.

   Returns:
       tuple: Variable length tuple with the arguments prepared to be able to be read by the bash interpreter.
   """

   result_list = []
   for arg in args:
      if type(arg) is list:
         if any(isinstance(x, (int, float)) for x in arg):
            list_with_strings = [str(x) for x in arg]
            joined_arg = ','.join(list_with_strings)
            result_list.append(joined_arg)
         else:
            joined_arg = ','.join(arg)
            result_list.append(joined_arg)
      elif type(arg) is dict:
         joined_arg = ','.join([f'{key}={value}' for key, value in arg.items()])
         result_list.append(joined_arg)
      else:
         result_list.append(str(arg))
   return tuple(result_list) if len(result_list) > 1 else result_list[0]

def msd_opt_prepper(nve_production_timestep_fs,
                    nve_production_duration_ns, 
                    reporting_frequency_coordinates_unwrapped,
                    fit_starting_percentage,
                    fit_ending_percentage) -> tuple:
   """Prepares the options for the MSD Analysis on the hpc.

   Args:
       nve_production_timestep_fs (int): The integration time step during the NVE production in femtoseconds.
       nve_production_duration_ns (int): The duration of the NVE production in nanoseconds.
       reporting_frequency_coordinates_unwrapped (int): The reporting frequency of the trajectory with unwrapped coordinates in simulation steps during the NVE production.
       fit_starting_percentage (int): Start of the linear fit in percent of NVE production duration.
       fit_ending_percentage (int): Ending of the linear fit in percent of NVE production duration.

   Returns:
       tuple: The options for the MSD Analysis on the hpc, fit_starting_frame, fit_ending_frame and time_between_frames.
   """
   
   time_between_frames = nve_production_timestep_fs * (10 ** (-3)) * reporting_frequency_coordinates_unwrapped
   conversion_factor_ps_to_frame = 1 / time_between_frames
   fit_starting_dec = fit_starting_percentage * (1 / 100)
   fit_ending_dec = fit_ending_percentage * (1 / 100)
   fit_starting_frame = int(fit_starting_dec * nve_production_duration_ns * (10 ** 3) * conversion_factor_ps_to_frame)
   fit_ending_frame = int(fit_ending_dec * nve_production_duration_ns * (10 ** 3) * conversion_factor_ps_to_frame)
   return fit_starting_frame, fit_ending_frame, time_between_frames