"""Main module.

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
import itertools
import json
import shutil
import re
import openmm.app.element as elem
from openmm.unit import *
from mixturemm.background import Simulationbox, Simulationsystem, Replica
import mixturemm.utils as m_utils


class Project:

    def __init__(self,
    workdir,
    simulation_platform='CUDA',
    simulation_properties={'DeviceIndex': '0', 'Precision': 'double'},
    total_number_molecules=[1000],
    init_box_side_length=[30],
    chi_water_s=[0,0.5,1],
    temperature_s=[298.15],
    npt_equilibration_pressure_s=[1],
    npt_equilibration_pressure_coupling_frequency=25,
    npt_equilibration_temperature_coupling_frequency=0.1,
    npt_equilibration_timestep_fs=2,
    npt_equilibration_duration_ns=10,
    reporting_frequency_state_npt_equilibration = 500,
    nvt_equilibration_temperature_coupling_frequency=0.1,
    nvt_equilibration_timestep_fs=2,
    nvt_equilibration_duration_ns=20,
    reporting_frequency_state_nvt_equilibration=2000,
    nve_production_timestep_fs=1,
    nve_production_duration_ns=20,
    reporting_frequency_coordinates_unwrapped=2000,
    reporting_frequency_coordinates_wrapped=4000,
    reporting_frequency_state_nve_production=4000,
    replica_count=10,
    pme_error_tolerance=0.000001,
    constraint_tolerance=0.0000001,
    cutoff_distance_nm=1.2,
    cutoff_switch_distance_nm=1.0
    ) -> None:
        """The project holds all parameters and functions to do a complete systematic simulation of different systems in openmm.

        Args:
            workdir (str): The path at which the project should be stored.
            simulation_platform (str, optional): The platform the simulation should run on. Defaults to 'CUDA'.
            simulation_properties (dict, optional): The platform properties for the simulation. Defaults to {'DeviceIndex': '0', 'Precision': 'double'}.
            total_number_molecules (list, optional): A list with the number of molecules that should be placed in the simulation box. Defaults to [1000].
            init_box_side_length (list, optional): The initial box side lengths that should be used by packmol. Defaults to [30].
            chi_water_s (list, optional): A list with all different molar fractions of water that should be simulated. Defaults to [0,0.5,1].
            temperature_s (list, optional): A list with all different temperatures in kelvin at which the boxes should be simulated. Defaults to [298.15].
            npt_equilibration_pressure_s (list, optional): A list with all pressures in bar at which the boxes should be simulated. Defaults to [1].
            npt_equilibration_pressure_coupling_frequency (int, optional): The pressure coupling frequency at which the Monte Carlo barostat should interact with the system and attempt a Monte Carlo move to adjust the volume during the NpT equilibration in simulation steps. Defaults to 25.
            npt_equilibration_temperature_coupling_frequency (float, optional): The friction coefficient of the Langevin thermostat during the NpT equilibration in inverse picoseconds. Defaults to 0.1.
            npt_equilibration_timestep_fs (int, optional): The integration time step during the NpT equilibration in femtoseconds. Defaults to 2.
            npt_equilibration_duration_ns (int, optional): The duration of the NpT equilibration in nanoseconds. Defaults to 10.
            reporting_frequency_state_npt_equilibration (int, optional): The reporting frequency of the openMM StateDataReporter in simulation steps during the NpT equilibration. Defaults to 500.
            nvt_equilibration_temperature_coupling_frequency (float, optional): The friction coefficient of the Langevin thermostat during the NVT equilibration in inverse picoseconds. Defaults to 0.1.
            nvt_equilibration_timestep_fs (int, optional): The integration time step during the NVT equilibration in femtoseconds. Defaults to 2.
            nvt_equilibration_duration_ns (int, optional): The duration of the NVT equilibration in nanoseconds. Defaults to 20.
            reporting_frequency_state_nvt_equilibration (int, optional): The reporting frequency of the openMM StateDataReporter in simulation steps during the NVT equilibration. Defaults to 2000.
            nve_production_timestep_fs (int, optional): The integration time step during the NVE production in femtoseconds. Defaults to 1.
            nve_production_duration_ns (int, optional): The duration of the NVE production in nanoseconds. Defaults to 20.
            reporting_frequency_coordinates_unwrapped (int, optional): The reporting frequency of the trajectory with unwrapped coordinates in simulation steps during the NVE production. Defaults to 2000.
            reporting_frequency_coordinates_wrapped (int, optional): The reporting frequency of the trajectory with wrapped coordinates in simulation steps during the NVE production. Defaults to 4000.
            reporting_frequency_state_nve_production (int, optional): The reporting frequency of the openMM StateDataReporter in simulation steps during the NVE production. Defaults to 4000.
            replica_count (int, optional): The number of times the simulation of one system should be replicated. Defaults to 10.
            pme_error_tolerance (float, optional): Decides how large the grid for PME is together with the cutoff. Defaults to 0.000001.
            constraint_tolerance (float, optional): The constraint tolerance specifies the distance tolerance within which constraints are maintained, as a fraction of the constrained distance. Defaults to 0.0000001
            cutoff_distance_nm (float, optional): The cutoff distance for nonbonded interactions in nanometers. Defaults to 1.2.
            cutoff_switch_distance_nm (float, optional): Starting point of the switching function that makes the energy go smoothly to 0 at the cutoff distance in nanometers. Defaults to 1.0.
        """
        
        self.workdir = m_utils.directory_maker(workdir)
        self.moldir = m_utils.directory_maker(f'{self.workdir}/molecules')
        self.boxdir = m_utils.directory_maker(f'{self.workdir}/boxes')
        self.resdir = m_utils.directory_maker(f'{self.workdir}/results')
        self.submitdir = m_utils.directory_maker(f'{self.workdir}/hpc_submission')
        self.forcedir = m_utils.directory_maker(f'{self.workdir}/forcefield')
        self.res_rawdir = m_utils.directory_maker(f'{self.resdir}/raw_data')
        self.res_partsdir = m_utils.directory_maker(f'{self.resdir}/parts')
        self.molecules = []
        self.simulation_boxes = []
        self.simulationsystems = []
        self.replicas = []
        self.molecule_names = []
        self.molecule_number_of_atoms = []
        self.molecule_abbreviations = []
        self.mixture_dict = {}
        self.forcefields = []
        self.elements = {}
        self.chi_water_s = chi_water_s
        self.total_number_molecules = total_number_molecules
        self.init_box_side_length = init_box_side_length
        self.water_name = 'Water'
        self.water_abbreviation = 'HOH'
        self.replica_number_dict = None
        self.simulation_platform = simulation_platform
        self.simulation_properties = simulation_properties
        self.temperature_s = temperature_s
        self.npt_equilibration_pressure_s = npt_equilibration_pressure_s
        self.npt_equilibration_pressure_coupling_frequency = npt_equilibration_pressure_coupling_frequency
        self.npt_equilibration_temperature_coupling_frequency = npt_equilibration_temperature_coupling_frequency
        self.npt_equilibration_timestep_fs = npt_equilibration_timestep_fs
        self.npt_equilibration_duration_ns = npt_equilibration_duration_ns
        self.nvt_equilibration_temperature_coupling_frequency = nvt_equilibration_temperature_coupling_frequency
        self.nvt_equilibration_timestep_fs = nvt_equilibration_timestep_fs
        self.nvt_equilibration_duration_ns = nvt_equilibration_duration_ns
        self.nve_production_timestep_fs = nve_production_timestep_fs
        self.nve_production_duration_ns = nve_production_duration_ns
        self.reporting_frequency_coordinates_unwrapped = reporting_frequency_coordinates_unwrapped
        self.reporting_frequency_coordinates_wrapped = reporting_frequency_coordinates_wrapped
        self.reporting_frequency_state_npt_equilibration = reporting_frequency_state_npt_equilibration
        self.reporting_frequency_state_nvt_equilibration = reporting_frequency_state_nvt_equilibration
        self.reporting_frequency_state_nve_production = reporting_frequency_state_nve_production
        self.replica_count = replica_count
        self.pme_error_tolerance = pme_error_tolerance
        self.constraint_tolerance = constraint_tolerance
        self.cutoff_distance_nm = cutoff_distance_nm
        self.cutoff_switch_distance_nm = cutoff_switch_distance_nm
        self.half_npt_equilibration_csv_columns = int((((self.npt_equilibration_duration_ns/self.npt_equilibration_timestep_fs)*(10**6))/reporting_frequency_state_npt_equilibration)/2)
        self.create_description()
        m_utils.add_chainsubmitter_script(self.submitdir)

    def create_description(self) -> None:
        """Creates a .json file with all input arguments of the project and info on all added files.
        """
            
        json_structure = m_utils.jsonize_project(self)
        if os.path.isfile(f'{self.workdir}/project_description.json') == False:
            with open(f'{self.workdir}/project_description.json', 'w', encoding='utf-8') as f:
                json.dump(json_structure, f, ensure_ascii=False, indent=4)
        else:
            m_utils.merge_project_descriptions(f'{self.workdir}/project_description.json', json_structure)
        self.description = f'{self.workdir}/project_description.json'

    def add_molecule(self, Molecule) -> None:
        """Adds a molecule to the project.

        Args:
            Molecule (object): Contains information on the molecule that is used for packing with packmol and to identify it.
        """

        if Molecule.use_as_water:
            self.water_name = Molecule.name
            self.water_abbreviation = Molecule.abbreviation
        else:
            self.molecules.append(Molecule)
            self.molecule_names.append(Molecule.name)
            self.molecule_number_of_atoms.append(Molecule.number_of_atoms)
            self.molecule_abbreviations.append(Molecule.abbreviation)
        json_structure = m_utils.jsonize_molecule(Molecule)
        m_utils.add_json_entry(f'{self.description}', json_structure, key='molecules', subkey=f'{Molecule.name}')

    def add_forcefield(self, Forcefield) -> None:
        """Adds a force field to the project.

        Args:
            Forcefield (object): Contains information on the force field.
        """

        self.forcefields.append(Forcefield.path)
        if Forcefield.elements is not None:
            self.elements.update(Forcefield.elements)
        json_structure = m_utils.jsonize_forcefield(Forcefield)
        m_utils.add_json_entry(f'{self.description}', json_structure, key='force fields', subkey=f'{Forcefield.ff_name}')

    def add_mixture(self, Mixture) -> None:
        """Adds a mixture to the project.

        Args:
            Mixture (object): Contains information on the mixture that is simulated.
        """

        self.mixture_dict = Mixture.mixture_dict

    def create_simulationboxes(self) -> None:
        """Creates and packs the simulation boxes.
        """

        LINKER = '-'

        for total_number_molecules, init_box_side_length in zip(self.total_number_molecules, self.init_box_side_length):
            for chi_water in self.chi_water_s:
                string_total_number_molecules = str(total_number_molecules)
                string_init_box_side_length = str(init_box_side_length)
                string_chi_water = str(chi_water)
                mixture_tuple = (string_total_number_molecules, string_init_box_side_length, string_chi_water)
                mixture_name = LINKER.join(mixture_tuple)
                simulation_box = Simulationbox(self, total_number_molecules, init_box_side_length, chi_water, self.water_name, mixture_name)
                self.simulation_boxes.append(simulation_box)
        for x in self.simulation_boxes:
            x.pack()
            x.conect_creator()

    def create_systems(self) -> None:
        """Assigns temperature, pressure and NpT equilibration parameters to the boxes.
        """

        parameter_list = [self.simulation_boxes, self.temperature_s, self.npt_equilibration_pressure_s]
        combination_list = list(itertools.product(*parameter_list))
        for combination in combination_list:
            simulationbox, temperature, pressure = combination
            system = Simulationsystem(self, simulationbox, temperature, pressure)
            self.simulationsystems.append(system)

    def create_replicas(self, start=1) -> None:
        """Creates replicas of the systems and assigns NVT equilibration and NVE production parameters to them.

        Args:
            start (int, optional): The number from which on the replica numbering will start. Defaults to 1.
        """

        self.replica_start = start
        replica_count = range(start, (self.replica_count + start))
        self.replica_number_dict = {}
        for tot in self.total_number_molecules:
            self.replica_number_dict[f'{tot}'] = self.replica_count
        parameter_list = [self.simulationsystems, replica_count]
        combination_list = list(itertools.product(*parameter_list))
        for combination in combination_list:
            system, replica_number = combination
            replica = Replica(self, system, replica_number)
            self.replicas.append(replica)

        m_utils.add_json_entry(f'{self.description}', start, key='project', subkey='simulation_parameters', subsubkey='replica_starting_number')

    def adjust_to_correct_density(self) -> None:
        """Runs the NpT equilibration to adjust the system volume and density to the correct value.
        """

        for system in self.simulationsystems:
            system.npt_equilibration()
    
    def simulate(self) -> None:
        """Runs the NVT equilibration and the NVE production.
        """

        for replica in self.replicas:
            replica.nvt_equilibration()
            replica.nve_production()

    def create_hpc_submission_adjust_to_correct_density(self,
    hpc_workspace='workspace',
    hpc_folder='workspace/project_name',
    hpc_scripts_folder='workspace/HPC',
    environment_name='openmm_new',
    scheduler='SBATCH',
    partition='accelerated',
    number_of_threads=152,
    number_of_gpus=4,
    chunk_size=4,
    max_number_of_jobs=None,
    max_runtime_hh_mm_ss='24:00:00',
    conda_module=False,
    chain_submission_number=1,
    dependency_type='afternotok',
    checkpoint_frequency=10000
    ) -> None:
        """Creates job submission scripts for bwUniCluster 2.0 or HoreKa containing the NpT equilibration.

        Args:
            hpc_workspace (str, optional): The path to the hpc workspace. Defaults to 'workspace'.
            hpc_folder (str, optional): The path to the project folder on the hpc. Defaults to 'workspace/project_name'.
            hpc_scripts_folder (str, optional): The path to the python scripts folder. Defaults to 'workspace/HPC'.
            environment_name (str, optional): The name of the anaconda environment containing the dependencies. Defaults to 'openmmfin'.
            scheduler (str, optional): The queue scheduler. Defaults to 'SBATCH'.
            partition (str, optional): The partition the job should be submitted to. Defaults to 'accelerated'.
            number_of_threads (int, optional): The number of requested threads. Defaults to 152.
            number_of_gpus (int, optional): The number of GPUs requested on a node. Defaults to 4.
            chunk_size (int, optional): Number of simulations to submit in a chunk, should correspond to the number of requested GPUs. Defaults to 4.
            max_number_of_jobs (int, optional): Maximum number of job submissions allowed on the hpc. Defaults to None.
            max_runtime_hh_mm_ss (str, optional): Maximum runtime for a job. Defaults to '24:00:00'.
            conda_module (bool): Adapts the header to whether there is a conda module available on the hpc or not. Defaults to False.
            chain_submission_number (int, optional): The number of job submissions that should be submitted as a sequential chain job. Defaults to 1.
            dependency_type (str, optional): Dependency keyword of the SLURM workload manager for chain jobs. Defaults to 'afternotok'.
            checkpoint_frequency (int, optional): The frequency at which checkpoints are saved to restart from once a job continues. Defaults to 10000.
        """

        forcefields_ = [re.sub(r'^.*?/forcefield', f'{hpc_folder}/forcefield', forcefield) for forcefield in self.forcefields]
        forcefields = ','.join(forcefields_)
        if self.elements:
            elements = m_utils.bash_arg_prepper(self.elements)
        else:
            elements = None
        if max_number_of_jobs is None:
            chunked_simulationsystems = [self.simulationsystems[i:i + chunk_size] for i in range(0, len(self.simulationsystems), chunk_size)]
        else:
            chunked_simulationsystems = [self.simulationsystems[i:i + chunk_size] for i in range(0, len(self.simulationsystems), chunk_size)][:max_number_of_jobs]
        for chunk in chunked_simulationsystems:
            GPU_INDEX_COUNTER = 0
            m_utils.hpc_submission_header(f'{self.submitdir}/chunk{id(chunk)}_adjust_to_correct_density.sh', 
                                          f'{id(chunk)}_adjust_to_correct_density', 
                                          hpc_workspace, 
                                          environment_name,
                                          scheduler,
                                          partition,
                                          number_of_threads,
                                          number_of_gpus, 
                                          max_runtime_hh_mm_ss,
                                          conda_module)
            for system in chunk:
                with open(f'{self.submitdir}/chunk{id(chunk)}_adjust_to_correct_density.sh', 'a', newline='\n') as file:
                    if 'DeviceIndex' and 'Precision' in self.simulation_properties:
                        simulation_properties_ = self.simulation_properties
                        simulation_properties_['DeviceIndex'] = f'{GPU_INDEX_COUNTER}'
                        simulation_properties = m_utils.bash_arg_prepper(simulation_properties_)
                    else:
                        simulation_properties_ = {'DeviceIndex': f'{GPU_INDEX_COUNTER}', 'Precision': 'double'}
                        simulation_properties = m_utils.bash_arg_prepper(simulation_properties_)
                    if self.simulation_platform == 'CUDA':
                        simulation_platform = self.simulation_platform
                    else:
                        simulation_platform = 'CUDA'
                    file.write(f'python {hpc_scripts_folder}/cluster_adjust_to_correct_density.py {hpc_folder} {simulation_properties} {simulation_platform} {forcefields} {system.pme_error_tolerance} {system.cutoff_distance_nm} {system.cutoff_switch_distance_nm} {system.name} {system.temperature} {system.pressure} {system.simulationbox_name} {system.npt_equilibration_pressure_coupling_frequency} {system.npt_equilibration_temperature_coupling_frequency} {system.npt_equilibration_timestep_fs} {system.npt_equilibration_steps} {system.reporting_frequency_state_npt_equilibration} {system.constraint_tolerance} {checkpoint_frequency} {elements} &\n')
                GPU_INDEX_COUNTER += 1
            with open(f'{self.submitdir}/chunk{id(chunk)}_adjust_to_correct_density.sh', 'a', newline='\n') as file:
                file.write('\nwait\nexit 0')
        with open(f'{self.submitdir}/submit_adjust_to_correct_density.sh', 'w', newline='\n') as file:
            file.write(f'''#!/bin/bash\n''')
            for chunk in chunked_simulationsystems:
                if chain_submission_number > 1:
                    file.write(f'''bash chainsubmitter.sh {chain_submission_number} {hpc_folder}/shfiles/chunk{id(chunk)}_adjust_to_correct_density.sh {dependency_type} {partition}\n''')
                else:
                    file.write(f'''sbatch chunk{id(chunk)}_adjust_to_correct_density.sh\n''')

    def create_hpc_submission_simulate(self,
    hpc_workspace='workspace',
    hpc_folder='workspace/project_name',
    hpc_scripts_folder='workspace/HPC',
    environment_name='openmm_new',
    scheduler='SBATCH',
    partition='accelerated',
    number_of_threads=152,
    number_of_gpus=4,
    chunk_size=4,
    max_number_of_jobs=None,
    max_runtime_hh_mm_ss='24:00:00',
    conda_module=False,
    chain_submission_number=1,
    dependency_type='afternotok',
    checkpoint_frequency=10000
    ) -> None:
        """Creates job submission scripts for bwUniCluster 2.0 containing the NVT equilibration and NVE production.

        Args:
            hpc_workspace (str, optional): The path to the hpc workspace. Defaults to 'workspace'.
            hpc_folder (str, optional): The path to the project folder on the hpc. Defaults to 'workspace/project_name'.
            hpc_scripts_folder (str, optional): The path to the python scripts folder. Defaults to 'workspace/HPC'.
            environment_name (str, optional): The name of the anaconda environment containing the dependencies. Defaults to 'openmmfin'.
            scheduler (str, optional): The queue scheduler. Defaults to 'SBATCH'.
            partition (str, optional): The partition the job should be submitted to. Defaults to 'accelerated'.
            number_of_threads (int, optional): The number of requested threads. Defaults to 152.
            number_of_gpus (int, optional): The number of GPUs requested on a node. Defaults to 4.
            chunk_size (int, optional): Number of simulations to submit in a chunk, should correspond to the number of requested GPUs. Defaults to 4.
            max_number_of_jobs (int, optional): Maximum number of job submissions allowed on the hpc. Defaults to None.
            max_runtime_hh_mm_ss (str, optional): Maximum runtime for a job. Defaults to '24:00:00'.
            conda_module (bool): Adapts the header to whether there is a conda module available on the hpc or not. Defaults to False.
            chain_submission_number (int, optional): The number of job submissions that should be submitted as a sequential chain job. Defaults to 1.
            dependency_type (str, optional): Dependency keyword of the SLURM workload manager for chain jobs. Defaults to 'afternotok'.
            checkpoint_frequency (int, optional): The frequency at which checkpoints are saved to restart from once a job continues. Defaults to 10000.
        """

        forcefields_ = [re.sub(r'^.*?/forcefield', f'{hpc_folder}/forcefield', forcefield) for forcefield in self.forcefields]
        forcefields = ','.join(forcefields_)
        if self.elements:
            elements = m_utils.bash_arg_prepper(self.elements)
        else:
            elements = None
        if max_number_of_jobs is None:
            chunked_replicas = [self.replicas[i:i + chunk_size] for i in range(0, len(self.replicas), chunk_size)]
        else:
            chunked_replicas = [self.replicas[i:i + chunk_size] for i in range(0, len(self.replicas), chunk_size)][:max_number_of_jobs]
        for chunk in chunked_replicas:
            GPU_INDEX_COUNTER = 0
            m_utils.hpc_submission_header(f'{self.submitdir}/chunk{id(chunk)}_simulate.sh', 
                                          f'{id(chunk)}_simulate', 
                                          hpc_workspace, 
                                          environment_name,
                                          scheduler,
                                          partition,
                                          number_of_threads,
                                          number_of_gpus,
                                          max_runtime_hh_mm_ss,
                                          conda_module)
            for replica in chunk:
                with open(f'{self.submitdir}/chunk{id(chunk)}_simulate.sh', 'a', newline='\n') as file:
                    if 'DeviceIndex' and 'Precision' in self.simulation_properties:
                        simulation_properties_ = self.simulation_properties
                        simulation_properties_['DeviceIndex'] = f'{GPU_INDEX_COUNTER}'
                        simulation_properties = m_utils.bash_arg_prepper(simulation_properties_)
                    else:
                        simulation_properties_ = {'DeviceIndex': f'{GPU_INDEX_COUNTER}', 'Precision': 'double'}
                        simulation_properties = m_utils.bash_arg_prepper(simulation_properties_)
                    if self.simulation_platform == 'CUDA':
                        simulation_platform = self.simulation_platform
                    else:
                        simulation_platform = 'CUDA'
                    file.write(f'python {hpc_scripts_folder}/cluster_simulate.py {hpc_folder} {simulation_properties} {simulation_platform} {forcefields} {replica.pme_error_tolerance} {replica.cutoff_distance_nm} {replica.cutoff_switch_distance_nm} {replica.name} {replica.replica_number} {replica.temperature} {replica.half_npt_equilibration_csv_columns} {replica.nvt_equilibration_temperature_coupling_frequency} {replica.nvt_equilibration_steps} {replica.nvt_equilibration_timestep_fs} {replica.reporting_frequency_state_nvt_equilibration} {replica.nve_production_timestep_fs} {replica.nve_production_steps} {replica.reporting_frequency_coordinates_unwrapped} {replica.reporting_frequency_coordinates_wrapped} {replica.reporting_frequency_state_nve_production} {checkpoint_frequency} {replica.constraint_tolerance} {elements} &\n')
                GPU_INDEX_COUNTER += 1
            with open(f'{self.submitdir}/chunk{id(chunk)}_simulate.sh', 'a', newline='\n') as file:
                file.write('\nwait\nexit 0')
        with open(f'{self.submitdir}/submit_simulate.sh', 'w', newline='\n') as file:
            file.write(f'''#!/bin/bash\n''')
            for chunk in chunked_replicas:
                if chain_submission_number > 1:
                    file.write(f'''bash chainsubmitter.sh {chain_submission_number} {hpc_folder}/shfiles/chunk{id(chunk)}_simulate.sh {dependency_type} {partition}\n''')
                else:
                    file.write(f'''sbatch chunk{id(chunk)}_simulate.sh\n''')

    def create_hpc_submission_analyze_density(self,
    hpc_workspace='workspace',
    hpc_folder='workspace/project_name',
    hpc_scripts_folder='workspace/HPC',
    environment_name='analysis',
    scheduler='SBATCH',
    partition='cpuonly',
    number_of_threads=152,
    max_runtime_hh_mm_ss='00:10:00',
    conda_module=False
    ) -> None:
        """Creates job submission scripts for bwUniCluster 2.0 containing the density anlysis.

        Args:
            hpc_workspace (str, optional): The path to the hpc workspace. Defaults to 'workspace'.
            hpc_folder (str, optional): The path to the project folder on the hpc. Defaults to 'workspace/project_name'.
            hpc_scripts_folder (str, optional): The path to the python scripts folder. Defaults to 'workspace/HPC'.
            environment_name (str, optional): The name of the anaconda environment containing the dependencies. Defaults to 'analysis'.
            scheduler (str, optional): The queue scheduler. Defaults to 'SBATCH'.
            partition (str, optional): The partition the job should be submitted to. Defaults to 'cpuonly'.
            number_of_threads (int, optional): The number of requested threads. Defaults to 152.
            max_runtime_hh_mm_ss (str, optional): Maximum runtime for a job. Defaults to '00:10:00'.
            conda_module (bool): Adapts the header to whether there is a conda module available on the hpc or not. Defaults to False.
        """

        description = re.sub(r'^.*?/project_description.json', f'{hpc_folder}/project_description.json', self.description)
        molecule_names, chi_water_s, temperature_s, npt_equilibration_pressure_s, name_list= m_utils.bash_arg_prepper(self.molecule_names, 
                                                                                                                      self.chi_water_s, 
                                                                                                                      self.temperature_s, 
                                                                                                                      self.npt_equilibration_pressure_s,
                                                                                                                      [system.name for system in self.simulationsystems])
        chunked_simulationsystems = [self.simulationsystems[i:i + 10] for i in range(0, len(self.simulationsystems), 10)]
        m_utils.hpc_submission_header(f'{self.submitdir}/id{id(self.simulationsystems)}_analyze_density.sh', 
                                      f'id{id(self.simulationsystems)}_analyze_density', 
                                      hpc_workspace, 
                                      environment_name,
                                      scheduler,
                                      partition,
                                      number_of_threads,
                                      None,
                                      max_runtime_hh_mm_ss,
                                      conda_module)
        for chunk in chunked_simulationsystems:
            for system in chunk:
                with open(f'{self.submitdir}/id{id(self.simulationsystems)}_analyze_density.sh', 'a', newline='\n') as file:
                    file.write(f'python {hpc_scripts_folder}/cluster_get_density.py {hpc_folder} {system.name} {self.half_npt_equilibration_csv_columns}  {system.total_number_molecules} {system.temperature} {system.pressure} {system.chi_water} {molecule_names} &\n')
            with open(f'{self.submitdir}/id{id(self.simulationsystems)}_analyze_density.sh', 'a', newline='\n') as file:
                file.write('wait\n')
        with open(f'{self.submitdir}/id{id(self.simulationsystems)}_analyze_density.sh', 'a', newline='\n') as file:
            file.write(f'python {hpc_scripts_folder}/cluster_join_density.py {chi_water_s} {temperature_s} {npt_equilibration_pressure_s} {hpc_folder} {name_list} {molecule_names} {self.water_name} {description} &\n')
            file.write('wait\nexit 0')
        with open(f'{self.submitdir}/submit_density_analysis.sh', 'w', newline='\n') as file:
            file.write(f'''#!/bin/bash\nsbatch id{id(self.simulationsystems)}_analyze_density.sh\n''')

    def create_hpc_submission_analyze_msd(self,
    hpc_workspace='workspace',
    hpc_folder='workspace/project_name',
    hpc_scripts_folder='workspace/HPC',
    environment_name='analysis',
    scheduler='SBATCH',
    partition='cpuonly',
    number_of_threads=152,
    parallel_running=5,
    submission_split=4,
    max_runtime_hh_mm_ss='24:00:00',
    conda_module=False,
    fit_starting_percentage=20,
    fit_ending_percentage=80,
    just_conclude=False
    ) -> None:
        """Creates job submission scripts for bwUniCluster 2.0 or HoreKa containing the mean squared displacement, self-diffusion coefficient and viscosity analysis.

        Args:
            hpc_workspace (str, optional): The path to the hpc workspace. Defaults to 'workspace'.
            hpc_folder (str, optional): The path to the project folder on the hpc. Defaults to 'workspace/project_name'.
            hpc_scripts_folder (str, optional): The path to the python scripts folder. Defaults to 'workspace/HPC'.
            environment_name (str, optional): The name of the anaconda environment containing the dependencies. Defaults to 'analysis'.
            scheduler (str, optional): The queue scheduler. Defaults to 'SBATCH'.
            partition (str, optional): The partition the job should be submitted to. Defaults to 'cpuonly'.
            number_of_threads (int, optional): The number of requested threads. Defaults to 152.
            parallel_running (int, optional): The number of replica analyses running in parallel. Defaults to 5.
            submission_split (int, optional): The number of job submissions the analysis is split in. Defaults to 4.
            max_runtime_hh_mm_ss (str, optional): Maximum runtime for a job. Defaults to '24:00:00'.
            conda_module (bool): Adapts the header to whether there is a conda module available on the hpc or not. Defaults to False.
            fit_starting_percentage (int, optional): Start of the linear fit in percent of NVE production duration. Defaults to 20.
            fit_ending_percentage (int, optional): Ending of the linear fit in percent of NVE production duration. Defaults to 80.
            just_conclude (bool, optional): If true, skips the generation of job submissions scripts and just generates the conclude script. Defaults to True.
        """

        description = re.sub(r'^.*?/project_description.json', f'{hpc_folder}/project_description.json', self.description)
        fit_starting_frame, fit_ending_frame, time_between_frames = m_utils.msd_opt_prepper(self.nve_production_timestep_fs,
                                                                                            self.nve_production_duration_ns,
                                                                                            self.reporting_frequency_coordinates_unwrapped,
                                                                                            fit_starting_percentage,
                                                                                            fit_ending_percentage)
        molecule_names, molecule_abbreviations, number_dict, total_number_molecules, chi_water_s, temperature_s, npt_equilibration_pressure_s, name_list = m_utils.bash_arg_prepper(self.molecule_names,
                                                                                                                                                                                    self.molecule_abbreviations,
                                                                                                                                                                                    self.replica_number_dict,
                                                                                                                                                                                    self.total_number_molecules,
                                                                                                                                                                                    self.chi_water_s,
                                                                                                                                                                                    self.temperature_s,
                                                                                                                                                                                    self.npt_equilibration_pressure_s,
                                                                                                                                                                                    [system.name for system in self.simulationsystems])
        chunked_replicas = [self.replicas[i:i + parallel_running] for i in range(0, len(self.replicas), parallel_running)]
        split_fac = int(len(chunked_replicas) / submission_split)
        submission_splits = [chunked_replicas[i:i + split_fac] for i in range(0, len(chunked_replicas), split_fac)]
        if len(submission_splits) > submission_split:
            split_fac = int(len(chunked_replicas) / submission_split) + 1
            submission_splits = [chunked_replicas[i:i + split_fac] for i in range(0, len(chunked_replicas), split_fac)]
        if not just_conclude:
            for split in submission_splits:
                m_utils.hpc_submission_header(f'{self.submitdir}/split{id(split)}_analyze_msd.sh', 
                                              f'split{id(split)}_analyze_msd', 
                                              hpc_workspace, 
                                              environment_name,
                                              scheduler,
                                              partition,
                                              number_of_threads,
                                              None,
                                              max_runtime_hh_mm_ss,
                                              conda_module)
                for chunk in split:
                    for replica in chunk:
                        with open(f'{self.submitdir}/split{id(split)}_analyze_msd.sh', 'a', newline='\n') as file:
                            file.write(f'python {hpc_scripts_folder}/cluster_get_sdc.py {hpc_folder} {replica.name} {replica.replica_number} {replica.chi_water} {self.water_abbreviation} {molecule_abbreviations} {time_between_frames} {fit_starting_frame} {fit_ending_frame} {replica.temperature} {replica.pressure} {replica.total_number_molecules} {molecule_names} {self.water_name} &\n')
                    with open(f'{self.submitdir}/split{id(split)}_analyze_msd.sh', 'a', newline='\n') as file:
                        file.write('wait\n')
                with open(f'{self.submitdir}/split{id(split)}_analyze_msd.sh', 'a', newline='\n') as file:
                    file.write('exit 0')
            with open(f'{self.submitdir}/submit_msd_analysis.sh', 'w', newline='\n') as file:
                file.write(f'''#!/bin/bash\n''')
                for split in submission_splits:
                    file.write(f'''sbatch split{id(split)}_analyze_msd.sh\n''')
        with open(f'{self.submitdir}/conclude_msd.sh', 'w', newline='\n') as file:
            file.write(f'''#!/bin/bash\n''')
            file.write(f'source {hpc_workspace}/conda/etc/profile.d/conda.sh\n')
            file.write(f'conda activate {environment_name}\n')
            file.write(f'python {hpc_scripts_folder}/cluster_join_sdc_and_comp_visc.py {total_number_molecules} {chi_water_s} {temperature_s} {npt_equilibration_pressure_s} {hpc_folder} {name_list} {molecule_names} {self.water_name} {number_dict} {molecule_abbreviations} {self.water_abbreviation} {description} {self.replica_start} &\n')
            file.write('wait\nexit 0')

    def create_hpc_submission_analyze_hbonds(self,
    hpc_workspace='workspace',
    hpc_folder='workspace/project_name',
    hpc_scripts_folder='workspace/HPC',
    environment_name='analysis',
    scheduler='SBATCH',
    partition='cpuonly',
    number_of_threads=152,
    parallel_running=5,
    submission_split=4,
    max_runtime_hh_mm_ss='24:00:00',
    conda_module=False,
    donors='donor_selection_string',
    hydrogens='hydrogen_selection_string',
    acceptors='acceptors_selection_string',
    just_conclude=False
    ) -> None:
        """Creates job submission scripts for bwUniCluster 2.0 or HoreKa containing the mean squared displacement, self-diffusion coefficient and viscosity analysis.

        Args:
            hpc_workspace (str, optional): The path to the hpc workspace. Defaults to 'workspace'.
            hpc_folder (str, optional): The path to the project folder on the hpc. Defaults to 'workspace/project_name'.
            hpc_scripts_folder (str, optional): The path to the python scripts folder. Defaults to 'workspace/HPC'.
            environment_name (str, optional): The name of the anaconda environment containing the dependencies. Defaults to 'analysis'.
            scheduler (str, optional): The queue scheduler. Defaults to 'SBATCH'.
            partition (str, optional): The partition the job should be submitted to. Defaults to 'cpuonly'.
            number_of_threads (int, optional): The number of requested threads. Defaults to 152.
            parallel_running (int, optional): The number of replica analyses running in parallel. Defaults to 5.
            submission_split (int, optional): The number of job submissions the analysis is split in. Defaults to 4.
            max_runtime_hh_mm_ss (str, optional): Maximum runtime for a job. Defaults to '24:00:00'.
            conda_module (bool): Adapts the header to whether there is a conda module available on the hpc or not. Defaults to False.
            donors (str, optional): Selection of all hydrogen bond donor atoms in MDAnalysis selection syntax. Defaults to 'donor_selection_string'.
            hydrogens (str, optional): Selection of all hydrogen bond hydrogen atoms in MDAnalysis selection syntax. Defaults to 'hydrogens_selection_string'.
            acceptors (str, optional): Selection of all hydrogen bond acceptor atoms in MDAnalysis selection syntax. Defaults to 'acceptors_selection_string'.
            just_conclude (bool, optional): If true, skips the generation of job submissions scripts and just generates the conclude script. Defaults to True.
        """
        
        description = re.sub(r'^.*?/project_description.json', f'{hpc_folder}/project_description.json', self.description)
        molecule_names, number_dict, chi_water_s, temperature_s, npt_equilibration_pressure_s, molecule_abbreviations, name_list = m_utils.bash_arg_prepper(self.molecule_names,
                                                                                                                                                                                                self.replica_number_dict,
                                                                                                                                                                                                self.chi_water_s,
                                                                                                                                                                                                self.temperature_s,
                                                                                                                                                                                                self.npt_equilibration_pressure_s,
                                                                                                                                                                                                self.molecule_abbreviations,
                                                                                                                                                                                                [system.name for system in self.simulationsystems])
        chunked_replicas = [self.replicas[i:i + parallel_running] for i in range(0, len(self.replicas), parallel_running)]
        split_fac = int(len(chunked_replicas) / submission_split)
        submission_splits = [chunked_replicas[i:i + split_fac] for i in range(0, len(chunked_replicas), split_fac)]
        if len(submission_splits) > submission_split:
            split_fac = int(len(chunked_replicas) / submission_split) + 1
            submission_splits = [chunked_replicas[i:i + split_fac] for i in range(0, len(chunked_replicas), split_fac)]
        if not just_conclude:
            for split in submission_splits:
                m_utils.hpc_submission_header(f'{self.submitdir}/split{id(split)}_analyze_hbonds.sh', 
                                              f'split{id(split)}_analyze_hbonds', 
                                              hpc_workspace, 
                                              environment_name,
                                              scheduler,
                                              partition,
                                              number_of_threads,
                                              None,
                                              max_runtime_hh_mm_ss,
                                              conda_module)
                for chunk in split:
                    for replica in chunk:
                        with open(f'{self.submitdir}/split{id(split)}_analyze_hbonds.sh', 'a', newline='\n') as file:
                            file.write(f'python {hpc_scripts_folder}/cluster_get_hbonds.py {hpc_folder} {replica.name} {replica.replica_number} {replica.chi_water} {replica.temperature} {replica.pressure} {replica.total_number_molecules} {molecule_names} "{donors}" "{hydrogens}" "{acceptors}" &\n')
                    with open(f'{self.submitdir}/split{id(split)}_analyze_hbonds.sh', 'a', newline='\n') as file:
                        file.write('wait\n')
                with open(f'{self.submitdir}/split{id(split)}_analyze_hbonds.sh', 'a', newline='\n') as file:
                    file.write('exit 0')
            with open(f'{self.submitdir}/submit_hbonds_analysis.sh', 'w', newline='\n') as file:
                file.write(f'''#!/bin/bash\n''')
                for split in submission_splits:
                    file.write(f'''sbatch split{id(split)}_analyze_hbonds.sh\n''')
        with open(f'{self.submitdir}/conclude_hbonds.sh', 'w', newline='\n') as file:
            file.write(f'''#!/bin/bash\n''')
            file.write(f'source {hpc_workspace}/conda/etc/profile.d/conda.sh\n')
            file.write(f'conda activate {environment_name}\n')
            file.write(f'python {hpc_scripts_folder}/cluster_join_hbonds.py {hpc_folder} {name_list} {number_dict} {description} {chi_water_s} {temperature_s} {npt_equilibration_pressure_s} {molecule_names} {self.water_name} {molecule_abbreviations} {self.water_abbreviation} {self.replica_start} &\n')
            file.write('wait\nexit 0')

    def hpc_extend_replica_folders(self, hpc_folder='workspace/project_name') -> None:
        """Creates newly added replica folders to the system folders on the hpc.

        Args:
            hpc_folder (str, optional): The path to the project folder on the hpc. Defaults to 'workspace/project_name'.
        """
        
        folder_list = []
        for replica in self.replicas:
            folder_list.append(re.sub(fr'^.*?/{replica.name}', f'{hpc_folder}/{replica.name}', replica.folder))
        with open(f'{self.submitdir}/extend_replica_folders.sh', 'w', newline='\n') as file:
            file.write(
f'''#!/bin/bash
''')
            for folder_string in folder_list:
                file.write(
f'''mkdir {folder_string}
''')

    def hpc_job_checker(self, 
    hpc_folder='workspace/project_name', 
    hpc_scripts_folder='workspace/HPC', 
    hpc_workspace='workspace', 
    environment_name='openmm_new', 
    verbose=True,
    energy_shift_tolerance_percent=1) -> None:
        """Creates a bash script to check the overall status of the simulations on the hpc. Checked properties are doneness percent, average temperature and total energy shift.

        Args:
            hpc_folder (str, optional): The path to the project folder on the hpc. Defaults to 'workspace/project_name'.
            hpc_scripts_folder (str, optional): The path to the python scripts folder. Defaults to 'workspace/HPC'.
            hpc_workspace (str, optional): The path to the hpc workspace. Defaults to 'workspace'.
            environment_name (str, optional): The name of the anaconda environment containing the dependencies. Defaults to 'openmm_new'.
            verbose (bool, optional): If true checks average temperature and total energy shift. Defaults to True.
            energy_shift_tolerance_percent (int, optional): The tolerated energy shift in percent. If the energy shift of a simulation is greater, a warning is issued. Defaults to 1.
        """

        system_file_list_ = []
        for system in self.simulationsystems:
            system_file_list_.append(re.sub(fr'^.*?/{system.name}', f'{hpc_folder}/{system.name}', system.folder))
        system_file_list = sorted(system_file_list_, key=lambda x: (float(x.split('/')[-1].split('-')[0]), float(x.split('/')[-1].split('-')[2]), float(x.split('/')[-1].split('-')[3]), float(x.split('/')[-1].split('-')[4])))
        with open(f'{self.submitdir}/system_list.txt', 'w', newline='\n') as file:
            for file_str in system_file_list:
                file.write(
f'''{file_str}
''')
        replica_file_list_ = []
        for replica in self.replicas:
            replica_file_list_.append(re.sub(fr'^.*?/{replica.name}', f'{hpc_folder}/{replica.name}', replica.folder))
        replica_file_list = sorted(replica_file_list_, key=lambda x: (float(x.split('/')[-2].split('-')[0]), float(x.split('/')[-2].split('-')[2]), float(x.split('/')[-2].split('-')[3]), float(x.split('/')[-2].split('-')[4]), float(x.split('/')[-1].split('_')[1])))
        replica_file_list.append('finalize')
        with open(f'{self.submitdir}/replica_list.txt', 'w', newline='\n') as file:
            for file_str in replica_file_list:
                file.write(
f'''{file_str}
''')
        with open(f'{self.submitdir}/job_checker.sh', 'w', newline='\n') as file:
            file.write(
f'''#!/bin/bash
source {hpc_workspace}/conda/etc/profile.d/conda.sh
conda activate {environment_name}
python {hpc_scripts_folder}/job_checker.py {self.nvt_equilibration_duration_ns} {self.nve_production_duration_ns} {verbose} {energy_shift_tolerance_percent} {self.npt_equilibration_duration_ns}
exit 0
''')

    def remove_done_jobs_from_systems(self, path_to_job_checks='workdir/job_checks.txt') -> None:
        """Removes all systems that are done according to the job checker from the project.

        Args:
            path_to_job_checks (str, optional): The path to the job checker output file. Defaults to 'workdir/job_checks.txt'.
        """
        
        with open(f'{path_to_job_checks}', 'r', newline='\n') as file:
            not_done_job_list_ = [line.rstrip() for line in file if 'not yet adjusted' in line]
            not_done_job_list = [item[:-17] for item in not_done_job_list_]
        filtered_systems = []
        for system in self.simulationsystems:
            for string in not_done_job_list:
                if system.name in string:
                    filtered_systems.append(system)
        self.simulationsystems = filtered_systems
    
    def remove_done_jobs_from_replicas(self, path_to_job_checks='workdir/job_checks.txt') -> None:
        """Removes all replicas that are done according to the job checker from the project.

        Args:
            path_to_job_checks (str, optional): The path to the job checker output file. Defaults to 'workdir/job_checks.txt'.
        """
        
        with open(f'{path_to_job_checks}', 'r', newline='\n') as file:
            not_done_job_list_ = [line.rstrip() for line in file if 'not yet finished' in line]
            not_done_job_list = [item[:-17] for item in not_done_job_list_]
        filtered_replicas = []
        for replica in self.replicas:
            for string in not_done_job_list:
                if replica.name in string and int(re.search(r'\d+$', string).group()) == replica.replica_number:
                    filtered_replicas.append(replica)
        self.replicas = filtered_replicas

    def overcharge_replicas(self, total_number_molecules=1000, overcharge_amount=10) -> None:
        """Adds the overcharge amount to the replicas of the systems specified by total number of molecules.

        Args:
            total_number_molecules (int, optional): Number of molecules that are placed in the simulation box. Defaults to 1000.
            overcharge_amount (int, optional): Amount of replicas that should be added. Defaults to 10.
        """
        
        filtered_systems = [system for system in self.simulationsystems if system.total_number_molecules in total_number_molecules]
        replica_count = range((self.replica_count + self.replica_start), (self.replica_count + self.replica_start + overcharge_amount))
        parameter_list = [filtered_systems, replica_count]
        combination_list = list(itertools.product(*parameter_list))
        for combination in combination_list:
            system, replica_number = combination
            replica = Replica(self, system, replica_number)
            self.replicas.append(replica)
        entry = f'with {overcharge_amount} replicas'
        m_utils.add_json_entry(f'{self.description}', entry, key='project', subkey='simulation_parameters', subsubkey='overcharged_boxes', subsubsubkey=f'{total_number_molecules} molecules')
        self.replica_number_dict[f'{total_number_molecules}'] = self.replica_count + self.replica_start + overcharge_amount

    def transfer_project(self,
    hpc_folder_old='workspace_old/project_name', 
    hpc_folder_new='workspace_new/project_name') -> None:
        """Transfers a project between two workspaces.

        Args:
            hpc_folder_old (str, optional): The path to the old project folder on the hpc. Defaults to 'workspace_old/project_name'.
            hpc_folder_new (str, optional): The path to the new project folder on the hpc. Defaults to 'workspace_new/project_name'.
        """

        system_file_list_old_ = []
        system_file_list_new_ = []
        for system in self.simulationsystems:
            system_file_list_old_.append(re.sub(fr'^.*?/{system.name}', f'{hpc_folder_old}/{system.name}', system.folder))
            system_file_list_new_.append(re.sub(fr'^.*?/{system.name}', f'{hpc_folder_new}/{system.name}', system.folder))
        system_file_list_old = sorted(system_file_list_old_, key=lambda x: (float(x.split('/')[-1].split('-')[0]), float(x.split('/')[-1].split('-')[2]), float(x.split('/')[-1].split('-')[3]), float(x.split('/')[-1].split('-')[4])))
        system_file_list_new = sorted(system_file_list_new_, key=lambda x: (float(x.split('/')[-1].split('-')[0]), float(x.split('/')[-1].split('-')[2]), float(x.split('/')[-1].split('-')[3]), float(x.split('/')[-1].split('-')[4])))
        with open(f'{self.submitdir}/transfer_project.py', 'w', newline='\n') as file:
            file.write('import shutil\n')
            for file_str_old, file_str_new in zip(system_file_list_old, system_file_list_new):
                file.write(
f'''shutil.copytree(r"{file_str_old}/npt_equilibration", r"{file_str_new}/npt_equilibration")
''')
        replica_file_list_old_ = []
        replica_file_list_new_ = []
        for replica in self.replicas:
            replica_file_list_old_.append(re.sub(fr'^.*?/{replica.name}', f'{hpc_folder_old}/{replica.name}', replica.folder))
            replica_file_list_new_.append(re.sub(fr'^.*?/{replica.name}', f'{hpc_folder_new}/{replica.name}', replica.folder))
        replica_file_list_old = sorted(replica_file_list_old_, key=lambda x: (float(x.split('/')[-2].split('-')[0]), float(x.split('/')[-2].split('-')[2]), float(x.split('/')[-2].split('-')[3]), float(x.split('/')[-2].split('-')[4]), float(x.split('/')[-1].split('_')[1])))
        replica_file_list_new = sorted(replica_file_list_new_, key=lambda x: (float(x.split('/')[-2].split('-')[0]), float(x.split('/')[-2].split('-')[2]), float(x.split('/')[-2].split('-')[3]), float(x.split('/')[-2].split('-')[4]), float(x.split('/')[-1].split('_')[1])))
        with open(f'{self.submitdir}/transfer_project.py', 'a', newline='\n') as file:
            for file_str_old, file_str_new in zip(replica_file_list_old, replica_file_list_new):
                file.write(
f'''shutil.copytree(r"{file_str_old}", r"{file_str_new}")
''')
        with open(f'{self.submitdir}/transfer_project.py', 'a', newline='\n') as file:
            file.write(
f'''shutil.copytree(r"{hpc_folder_old}/boxes", r"{hpc_folder_new}/boxes")
shutil.copytree(r"{hpc_folder_old}/forcefield", r"{hpc_folder_new}/forcefield")
shutil.copytree(r"{hpc_folder_old}/hpc_submission", r"{hpc_folder_new}/hpc_submission")
shutil.copytree(r"{hpc_folder_old}/molecules", r"{hpc_folder_new}/molecules")
shutil.copytree(r"{hpc_folder_old}/results", r"{hpc_folder_new}/results")
''')


class Molecule:

    def __init__(self,
    workdir,
    name,
    number_of_atoms,
    path_to_pdb,
    abbreviation,
    smiles,
    inchi,
    molar_mass,
    use_as_water=False
    ) -> None:
        """Contains information on the molecule that is used for packing with packmol and to identify it.

        Args:
            workdir (str): The path at which the project is stored.
            name (str): Name of the molecule.
            number_of_atoms (int): The number of atoms that the molecule consists of.
            path_to_pdb (str): The path to the molecule PDB file.
            abbreviation (str): The abbreviation of the molecule that is used inside the PDB file.
            smiles (str): The SMILES code of the molecule.
            inchi (str): The InChI code of the molecule.
            molar_mass (str): The molar mass of the molecule.
            use_as_water (bool, optional): Decides whether the molecule is used as water or not. If used as water, the water molar fraction is applied to it while packing boxes. Defaults to False.
        """
        
        self.moldir = f'{workdir}/molecules'
        self.name = name
        self.number_of_atoms = number_of_atoms
        self.abbreviation = abbreviation
        self.smiles = smiles
        self.inchi = inchi
        self.molar_mass = molar_mass
        self.use_as_water = False
        if use_as_water:
            self.use_as_water = True
        shutil.copy(path_to_pdb, f'{self.moldir}/{self.name}.pdb')
        self.path = f'{self.moldir}/{self.name}.pdb'


class Mixture:

    def __init__(self,
    mixture_dict
    ) -> None:
        """Contains information on the mixture that is simulated.

        Args:
            mixture_dict (dict): Dictionary with molecule names and their corresponding mole fraction in the mixture.
        """
        
        self.mixture_dict = mixture_dict


class Forcefield:

    def __init__(self, 
    workdir,
    path_to_XML, 
    ff_name, 
    built_in=True, 
    elements=None) -> None:
        """Contains information on the force field.

        Args:
            workdir (str): The path at which the project is stored.
            path_to_XML (str): The path to the force field XML file.
            ff_name (str): The name under which the force field should be stored.
            built_in (bool, optional): Decides whether the fore field is part of openMM or from an external source. If true, declares that a built in forcefield of openmm will be used. Defaults to True.
            elements (dict, optional): OpenMM uses Elements to match molecules to the force field. If atomtypes in the force field are unknown to OpenMM, they are registered by giving their name and mass in the form {'name' : mass}. Defaults to None.
        """
        
        self.forcedir = f'{workdir}/forcefield'
        self.ff_name = ff_name
        self.built_in = built_in
        self.elements = elements
        if built_in == False:
            shutil.copy(path_to_XML, f'{self.forcedir}/{ff_name}.xml')
            self.path = f'{self.forcedir}/{ff_name}.xml'
        else:
            self.path = path_to_XML
        if elements is not None:
            for name, mass in elements.items():
                _ = elem.Element(number=0, name=name, symbol=name, mass=mass*amu)