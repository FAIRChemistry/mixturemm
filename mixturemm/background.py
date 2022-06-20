"""Module with indirectly used classes.

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

import itertools
import os
import subprocess
from openmm.app import *
from openmm import *
from openmm.unit import *
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from mixturemm.utils import directory_maker


class Simulationbox:

    def __init__(self,
    Project,
    total_number_molecules,
    init_box_side_length,
    chi_water,
    water_name,
    mixture_name
    ):
        """The simulationbox contains all information that is needed to pack each simulation box in packmol.

        Args:
            Project (object): The project holds all parameters and functions to do a complete systematic simulation of different systems in openmm.
            total_number_molecules (int): A list with the number of molecules that should be placed in the simulation box.
            init_box_side_length (int): The initial box side length that should be used by packmol.
            chi_water (flot): The molar fraction of water of this simulation box.
            water_name (str): The name that was assigned to the water molecule.
            mixture_name (str): The name created from a tuple of total_number_molecules, init_box_side_length and chi_water.
        """
        
        self.boxdir = Project.boxdir
        self.moldir = Project.moldir
        self.molecule_names = Project.molecule_names
        self.molecule_number_of_atoms = Project.molecule_number_of_atoms
        self.mixture_dict = Project.mixture_dict
        self.total_number_molecules = total_number_molecules
        self.init_box_side_length = init_box_side_length
        self.chi_water = chi_water
        self.water_name = water_name
        self.mixture_name = mixture_name
        self.path = f'{self.boxdir}/mixture_{mixture_name}.pdb'
        self.skip_bonds = False

    def pack(self):
        """Creates an input file for packmol with the initial box side length and the exact number of each involved molecule and runs packmol with it.
        """

        if os.path.isfile(f'{self.path}'):
            self.skip_bonds = True
        else:
            mixture_name = self.mixture_name
            operating_folder = self.boxdir
            moldir = self.moldir
            molecule_number_of_atoms = self.molecule_number_of_atoms
            total_number_molecules = int(self.total_number_molecules)
            init_box_side_length = int(self.init_box_side_length)
            mixture_dict = self.mixture_dict
            molecule_type_beginning_count = [0]
            COUNTER = 0
            chi_water = float(self.chi_water)
            water_number_molecules_ = chi_water * total_number_molecules
            water_number_molecules = int(water_number_molecules_)
            remaining_number_molecules = total_number_molecules - water_number_molecules
            with open(f'{operating_folder}/mixture_{mixture_name}.inp', 'w') as file:
                file.write(f'tolerance 2.0\nfiletype pdb\noutput {operating_folder}/mixture_{mixture_name}.pdb\n\n')
            with open(f'{operating_folder}/mixture_{mixture_name}.inp', 'a') as file:
                if chi_water < 1:
                    for molecule_name in mixture_dict:
                        number = int(round(mixture_dict[molecule_name] * remaining_number_molecules))
                        file.write(f'structure {moldir}/{molecule_name}.pdb\n\tnumber {number}\n\tinside box 0. 0. 0. {init_box_side_length}. {init_box_side_length}. {init_box_side_length}.\nend structure\n\n')
                if water_number_molecules == 0:
                    number_each_molecule_ = [mixture_dict[molecule_name]*remaining_number_molecules for molecule_name in mixture_dict]
                    number_each_molecule = [int(x) for x in number_each_molecule_]
                    for x, y in zip(number_each_molecule, molecule_number_of_atoms):
                        molecule_type_beginning_count.append(sum(molecule_type_beginning_count[COUNTER:]) + x * y)
                        COUNTER += 1
                    del molecule_type_beginning_count[-1]
                    self.molecule_type_beginning_count = molecule_type_beginning_count
                    self.number_each_molecule = number_each_molecule
                else:
                    file.write(f'structure {moldir}/{self.water_name}.pdb\n\tnumber {water_number_molecules}\n\tinside box 0. 0. 0. {init_box_side_length}. {init_box_side_length}. {init_box_side_length}.\nend structure\n\n')
                    number_each_molecule_ = [mixture_dict[molecule_name]*remaining_number_molecules for molecule_name in mixture_dict] + [water_number_molecules]
                    number_each_molecule = [int(x) for x in number_each_molecule_]
                    for x, y in zip(number_each_molecule, molecule_number_of_atoms):
                        molecule_type_beginning_count.append(sum(molecule_type_beginning_count[COUNTER:]) + x * y)
                        COUNTER += 1
                    self.molecule_type_beginning_count = molecule_type_beginning_count
                    self.number_each_molecule = number_each_molecule
            with open(f'{operating_folder}/mixture_{mixture_name}.inp', 'a') as file:
                file.write(f'add_box_sides 1.0')
            arguments = f'packmol < {operating_folder}/mixture_{mixture_name}.inp'
            _ = subprocess.run(arguments, shell= True)

    def conect_creator(self):
        """Takes the CONECT entries of all involved molecules from their PDB files and creates the CONECT entry for the simulation box accordingly.
        """

        if self.skip_bonds:
            pass
        else:
            path = self.path
            moldir = self.moldir
            molecule_names = self.molecule_names
            number_each_molecule = self.number_each_molecule
            number_each_molecule_dict = dict(zip(molecule_names, number_each_molecule))
            molecule_atnums = self.molecule_number_of_atoms
            molecule_atnums_dict = dict(zip(molecule_names, molecule_atnums))
            molecule_type_beginning_count = self.molecule_type_beginning_count
            molecule_type_beginning_count_dict = dict(zip(molecule_names, molecule_type_beginning_count))
            COUNTER = 0
            conect_raw_dict = dict.fromkeys(molecule_names)
            for name, _ in conect_raw_dict.items():
                conect_raw_ = []
                for line in open(f'{moldir}/{name}.pdb'):
                    if line[:6] == 'CONECT':
                        conect_raw_.append(line)
                conect_raw_string = ''.join(conect_raw_)
                conect_raw = [int(s) for s in conect_raw_string.split() if s.isdigit()]
                conect_raw_dict[name] = conect_raw
            conect_extended_dict = dict.fromkeys(molecule_names)
            for name, _ in conect_extended_dict.items():
                conect_extend = []
                range_index = number_each_molecule_dict[name]
                for _ in range(range_index):
                    conect_extend.append([])
                for x in range(range_index):
                    conect_extend[x].append([(f+((molecule_atnums_dict[name])*x)+molecule_type_beginning_count_dict[name]) for f in conect_raw_dict[name]])
                conect_extended_dict[name] = conect_extend
            conect_frame_dict = dict.fromkeys(molecule_names)
            for name, _ in conect_frame_dict.items():
                conect_frame = []
                for line in open(f'{moldir}/{name}.pdb'):
                    column_1 = line[12:16].strip()
                    column_2 = line[17:21].strip()
                    column_3 = line[22:26].strip()
                    column_4 = line[27:31].strip()
                    if line[:6] == 'CONECT' and column_4.isdigit() == True:
                        conect_frame.append(5)
                    elif line[:6] == 'CONECT' and column_3.isdigit() == True:
                        conect_frame.append(4)
                    elif line[:6] == 'CONECT' and column_2.isdigit() == True:
                        conect_frame.append(3)
                    elif line[:6] == 'CONECT' and column_1.isdigit() == True:
                        conect_frame.append(2)
                    else:
                        pass
                conect_frame_dict[name] = conect_frame
            conect_extended_flat_dict = dict.fromkeys(molecule_names)
            for name, _ in conect_extended_flat_dict.items():
                flattend_ = list(itertools.chain.from_iterable(conect_extended_dict[name]))
                flattend = list(itertools.chain.from_iterable(flattend_))
                conect_extended_flat_dict[name] = flattend
            conect_frame_extended_dict = dict.fromkeys(molecule_names)
            for name, _ in conect_frame_extended_dict.items():
                conect_frame_extend = []
                for __ in itertools.repeat(None, number_each_molecule_dict[name]):
                    conect_frame_extend.extend(conect_frame_dict[name])
                conect_frame_extended_dict[name] = conect_frame_extend
            conect_cleaned = {k:v for k, v in conect_extended_flat_dict.items() if v is not None}
            conect_frame_cleaned = {k:v for k, v in conect_frame_extended_dict.items() if v is not None}
            conect_final = []
            conect_frame_final = []
            for name, _ in conect_cleaned.items():
                conect_final.extend(conect_cleaned[name])
            for name, _ in conect_frame_cleaned.items():
                conect_frame_final.extend(conect_frame_cleaned[name])
            file = open(f'{path}', 'rt')
            data = file.read()
            data = data.replace('END\n', '')
            file.close()
            file = open(f'{path}', 'wt')
            file.write(data)
            file.close()
            with open (f'{path}', 'a') as file:
                for identifier, _ in itertools.zip_longest(conect_frame_final, conect_final):
                    if identifier == 1:
                        file.write(f'CONECT{conect_final[COUNTER]:5}\n')
                        COUNTER += 1
                    elif identifier == 2:
                        file.write(f'CONECT{conect_final[COUNTER]:5}{conect_final[COUNTER + 1]:5}\n')
                        COUNTER += 2
                    elif identifier == 3:
                        file.write(f'CONECT{conect_final[COUNTER]:5}{conect_final[COUNTER + 1]:5}{conect_final[COUNTER + 2]:5}\n')
                        COUNTER += 3
                    elif identifier == 4:
                        file.write(f'CONECT{conect_final[COUNTER]:5}{conect_final[COUNTER + 1]:5}{conect_final[COUNTER + 2]:5}{conect_final[COUNTER + 3]:5}\n')
                        COUNTER += 4
                    elif identifier == 5:
                        file.write(f'CONECT{conect_final[COUNTER]:5}{conect_final[COUNTER + 1]:5}{conect_final[COUNTER + 2]:5}{conect_final[COUNTER + 3]:5}{conect_final[COUNTER + 4]:5}\n')
                        COUNTER += 5
                    else:
                        pass
                file.write('END')


class Simulationsystem:

    def __init__(self,
    Project,
    Simulationbox,
    temperature,
    pressure,
    ):
        """The simulationsystem adds the respective temperature and pressure to each box and contains the parameters for the NpT equilibration that is done to adjust the system to the correct density.

        Args:
            Project (object): The project holds all parameters and functions to do a complete systematic simulation of different systems in openmm.
            Simulationbox (object): The simulationbox contains all information that is needed to pack each simulation box in packmol.
            temperature (float): Temperature in kelvin at which the system should be simulated.
            pressure (float): Pressures in bar at which the system should be simulated.
        """
        LINKER = '-'
        self.temperature = temperature
        self.pressure = pressure
        self.total_number_molecules = Simulationbox.total_number_molecules
        self.simulationbox_name = Simulationbox.mixture_name
        temperature_string = str(temperature)
        pressure_string = str(pressure)
        combined_strings = (self.simulationbox_name, temperature_string, pressure_string)
        system_name = LINKER.join(combined_strings)
        self.name = f'{system_name}'
        self.folder = f'{Project.workdir}/{system_name}'
        self.npt_equilibration_subfolder = f'{self.folder}/npt_equilibration'
        self.simulationbox_path = Simulationbox.path
        self.chi_water = Simulationbox.chi_water
        directory_maker(self.folder)
        directory_maker(self.npt_equilibration_subfolder)
        self.npt_equilibration_pressure_s = Project.npt_equilibration_pressure_s
        self.npt_equilibration_pressure_coupling_frequency = Project.npt_equilibration_pressure_coupling_frequency
        self.npt_equilibration_temperature_coupling_frequency = Project.npt_equilibration_temperature_coupling_frequency
        self.npt_equilibration_timestep_fs = Project.npt_equilibration_timestep_fs
        self.npt_equilibration_duration_ns = Project.npt_equilibration_duration_ns
        self.reporting_frequency_state_npt_equilibration = Project.reporting_frequency_state_npt_equilibration
        self.cutoff_distance_nm = Project.cutoff_distance_nm
        self.cutoff_switch_distance_nm = Project.cutoff_switch_distance_nm
        self.simulation_platform = Project.simulation_platform
        self.forcefields = Project.forcefields
        self.npt_equilibration_steps = int((self.npt_equilibration_duration_ns/self.npt_equilibration_timestep_fs) * (10 ** 6))
        self.openmm_forcefield = ForceField(*self.forcefields)
        self.openmm_platform = Platform.getPlatformByName(f'{self.simulation_platform}')
        self.openmm_properties = Project.simulation_properties
        self.pme_error_tolerance = Project.pme_error_tolerance
        self.constraint_tolerance = Project.constraint_tolerance
    
    def npt_equilibration(self):
        """Runs a NpT equilibration to adjust the box to the correct density using a Monte Carlo barostat.
        """

        START_TEMPERATURE = 1
        TEMPERATURE_STEP = 0.1
        heating_interval_steps = int(1000 / self.npt_equilibration_timestep_fs)
        folder = self.npt_equilibration_subfolder
        openmm_properties = self.openmm_properties
        openmm_platform = self.openmm_platform
        pme_error_tolerance = self.pme_error_tolerance
        cutoff_distance_nm = self.cutoff_distance_nm
        cutoff_switch_distance_nm = self.cutoff_switch_distance_nm
        temperature = self.temperature
        pressure = self.pressure
        openmm_forcefield = self.openmm_forcefield
        simulationbox_path = self.simulationbox_path
        npt_equilibration_pressure_coupling_frequency = self.npt_equilibration_pressure_coupling_frequency
        npt_equilibration_temperature_coupling_frequency = self.npt_equilibration_temperature_coupling_frequency
        npt_equilibration_timestep_ps = self.npt_equilibration_timestep_fs * (10 ** (-3))
        npt_equilibration_steps = self.npt_equilibration_steps
        reporting_frequency_state_npt_equilibration = self.reporting_frequency_state_npt_equilibration
        constraint_tolerance = self.constraint_tolerance

        topology_npt_equilibration = PDBFile(f'{simulationbox_path}')
        system_npt_equilibration = openmm_forcefield.createSystem(topology_npt_equilibration.topology, 
                                                                  nonbondedMethod=PME, 
                                                                  nonbondedCutoff=cutoff_distance_nm * nanometer,
                                                                  switchDistance=cutoff_switch_distance_nm * nanometer,
                                                                  constraints=HBonds,
                                                                  rigidWater=True,
                                                                  ewaldErrorTolerance=pme_error_tolerance,
                                                                  useDispersionCorrection=True)
        integrator_npt_equilibration = LangevinMiddleIntegrator(START_TEMPERATURE * kelvin, 10 / picosecond, npt_equilibration_timestep_ps * picoseconds)
        integrator_npt_equilibration.setConstraintTolerance(constraint_tolerance)
        barostat_monte_carlo = MonteCarloBarostat(pressure * bar, temperature * kelvin, 0)
        system_npt_equilibration.addForce(barostat_monte_carlo)
        simulation_npt_equilibration = Simulation(topology_npt_equilibration.topology, 
                                                  system_npt_equilibration, 
                                                  integrator_npt_equilibration, 
                                                  openmm_platform, 
                                                  openmm_properties)
        simulation_npt_equilibration.context.setPositions(topology_npt_equilibration.positions)
        simulation_npt_equilibration.minimizeEnergy(tolerance=0.1 * kilojoule / mole, maxIterations=50000)
        simulation_npt_equilibration.context.setVelocitiesToTemperature(START_TEMPERATURE * kelvin)
        simulation_npt_equilibration.reporters.append(StateDataReporter(f'{folder}/npt_equilibration.csv',
                                                                        reporting_frequency_state_npt_equilibration,
                                                                        time=True, 
                                                                        totalEnergy=True, 
                                                                        density=True, 
                                                                        volume=True,  
                                                                        potentialEnergy=True, 
                                                                        kineticEnergy=True, 
                                                                        temperature=True))
        for temp in np.arange(START_TEMPERATURE, temperature, TEMPERATURE_STEP):
            integrator_npt_equilibration.setTemperature(temp * kelvin)
            simulation_npt_equilibration.step(heating_interval_steps)
        simulation_npt_equilibration.step(heating_interval_steps * 10)
        integrator_npt_equilibration.setTemperature(temperature * kelvin)
        integrator_npt_equilibration.setFriction(npt_equilibration_temperature_coupling_frequency)
        barostat_monte_carlo.setFrequency(npt_equilibration_pressure_coupling_frequency)
        simulation_npt_equilibration.step(npt_equilibration_steps)
        system_npt_equilibration.removeForce(system_npt_equilibration.getNumForces() - 1)
        simulation_npt_equilibration.saveState(f'{folder}/npt_equilibration_state.xml')
        endpositions_npt_equilibration = simulation_npt_equilibration.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(topology_npt_equilibration.topology, endpositions_npt_equilibration, open(f'{folder}/npt_equilibration_end.pdb', 'w'))
        df = pd.read_csv(f'{folder}/npt_equilibration.csv')
        df.columns = ["Time (ps)","Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)","Temperature (K)","Box Volume (nm^3)","Density (g/mL)"]
        fig, axes = plt.subplots(2, 3, figsize=(40, 20))
        sns.scatterplot(data=df, ax=axes[0,0], x="Time (ps)", y="Box Volume (nm^3)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[0,1], x="Time (ps)", y="Density (g/mL)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[0,2], x="Time (ps)", y="Total Energy (kJ/mole)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[1,0], x="Time (ps)", y="Temperature (K)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[1,1], x="Time (ps)", y="Potential Energy (kJ/mole)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[1,2], x="Time (ps)", y="Kinetic Energy (kJ/mole)", s=2, marker='o')
        plt.savefig(f'{folder}/npt_equilibration.svg', dpi=300)
        plt.savefig(f'{folder}/npt_equilibration.png', dpi=300)


class Replica:

    def __init__(self,
    Project,
    Simulationsystem,
    replica_number
    ):
        """The replica of a simulationsystem contains the parameters for the NVT equilibration and the NVE production.

        Args:
            Project (object): The project holds all parameters and functions to do a complete systematic simulation of different systems in openmm.
            Simulationsystem (object): The simulationsystem adds the respective temperature and pressure to each box and contains the parameters for the NpT equilibration that is done to adjust the system to the correct density.
            replica_number (int): The number assigned to this replica.
        """

        self.folder = f'{Simulationsystem.folder}/replica_{replica_number}'
        directory_maker(self.folder)
        self.npt_equilibration_subfolder = Simulationsystem.npt_equilibration_subfolder
        self.name = Simulationsystem.name
        self.simulationbox_path = Simulationsystem.simulationbox_path
        self.replica_number = replica_number
        self.chi_water = Simulationsystem.chi_water
        self.temperature = Simulationsystem.temperature
        self.pressure = Simulationsystem.pressure
        self.total_number_molecules = Simulationsystem.total_number_molecules
        self.simulationbox_name = Simulationsystem.simulationbox_name
        self.nvt_equilibration_temperature_coupling_frequency = Project.nvt_equilibration_temperature_coupling_frequency
        self.nvt_equilibration_timestep_fs = Project.nvt_equilibration_timestep_fs
        self.nvt_equilibration_duration_ns = Project.nvt_equilibration_duration_ns
        self.nve_production_timestep_fs = Project.nve_production_timestep_fs
        self.nve_production_duration_ns = Project.nve_production_duration_ns
        self.reporting_frequency_coordinates_unwrapped = Project.reporting_frequency_coordinates_unwrapped
        self.reporting_frequency_coordinates_wrapped = Project.reporting_frequency_coordinates_wrapped
        self.reporting_frequency_state_nvt_equilibration = Project.reporting_frequency_state_nvt_equilibration
        self.reporting_frequency_state_nve_production = Project.reporting_frequency_state_nve_production
        self.cutoff_distance_nm = Project.cutoff_distance_nm
        self.cutoff_switch_distance_nm = Project.cutoff_switch_distance_nm
        self.simulation_platform = Project.simulation_platform
        self.half_npt_equilibration_csv_columns = Project.half_npt_equilibration_csv_columns
        self.forcefields = Project.forcefields
        self.nvt_equilibration_steps = int((self.nvt_equilibration_duration_ns/self.nvt_equilibration_timestep_fs) * (10 ** 6))
        self.nve_production_steps = int((self.nve_production_duration_ns/self.nve_production_timestep_fs) * (10 ** 6))
        self.openmm_forcefield = ForceField(*self.forcefields)
        self.openmm_platform = Platform.getPlatformByName(f'{self.simulation_platform}')
        self.openmm_properties = Project.simulation_properties
        self.pme_error_tolerance = Project.pme_error_tolerance
        self.constraint_tolerance = Project.constraint_tolerance

    def nvt_equilibration(self):
        """Runs a NVT equilibration of the system.
        """

        folder = self.folder
        npt_equilibration_subfolder = self.npt_equilibration_subfolder
        npt_equilibration_state = f'{npt_equilibration_subfolder}/npt_equilibration_state.xml'
        npt_equilibration_end_simulationbox = f'{npt_equilibration_subfolder}/npt_equilibration_end.pdb'
        openmm_properties = self.openmm_properties
        openmm_platform = self.openmm_platform
        pme_error_tolerance = self.pme_error_tolerance
        cutoff_distance_nm = self.cutoff_distance_nm
        cutoff_switch_distance_nm = self.cutoff_switch_distance_nm
        temperature = self.temperature
        openmm_forcefield = self.openmm_forcefield
        nvt_equilibration_temperature_coupling_frequency = self.nvt_equilibration_temperature_coupling_frequency
        nvt_equilibration_timestep_ps = self.nvt_equilibration_timestep_fs * (10 ** (-3))
        nvt_equilibration_steps = self.nvt_equilibration_steps
        half_npt_equilibration_csv_columns = self.half_npt_equilibration_csv_columns
        reporting_frequency_state_nvt_equilibration = self.reporting_frequency_state_nvt_equilibration
        constraint_tolerance = self.constraint_tolerance

        topology_nvt_equilibration = PDBFile(f'{npt_equilibration_end_simulationbox}')      
        system_nvt_equilibration = openmm_forcefield.createSystem(topology_nvt_equilibration.topology, 
                                                                  nonbondedMethod=PME, 
                                                                  nonbondedCutoff=cutoff_distance_nm * nanometer,
                                                                  switchDistance=cutoff_switch_distance_nm * nanometer,
                                                                  constraints=HBonds,
                                                                  rigidWater=True,
                                                                  ewaldErrorTolerance=pme_error_tolerance,
                                                                  useDispersionCorrection=True)
        integrator_nvt_equilibration = LangevinMiddleIntegrator(temperature * kelvin, 
                                                            nvt_equilibration_temperature_coupling_frequency / picosecond, 
                                                            nvt_equilibration_timestep_ps * picoseconds)
        integrator_nvt_equilibration.setConstraintTolerance(constraint_tolerance)
        simulation_nvt_equilibration = Simulation(topology_nvt_equilibration.topology, 
                                                  system_nvt_equilibration, 
                                                  integrator_nvt_equilibration, 
                                                  openmm_platform, 
                                                  openmm_properties,
                                                  state=npt_equilibration_state)
        npt_equilibration_df = pd.read_csv(f'{npt_equilibration_subfolder}/npt_equilibration.csv')
        half_npt_equilibration_df = npt_equilibration_df.tail(n=half_npt_equilibration_csv_columns)
        volume_npt_equilibration = half_npt_equilibration_df.loc[:, 'Box Volume (nm^3)']
        average_volume_npt_equilibration = volume_npt_equilibration.mean()
        estimated_new_side_length = (average_volume_npt_equilibration ** (1 / 3))
        positions = simulation_nvt_equilibration.context.getState(getPositions=True).getPositions()
        simulation_nvt_equilibration.context.reinitialize()
        simulation_nvt_equilibration.context.setPositions(positions)
        simulation_nvt_equilibration.context.setPeriodicBoxVectors(Vec3(x=estimated_new_side_length, y=0.0, z=0.0), 
                                                                   Vec3(x=0.0, y=estimated_new_side_length, z=0.0), 
                                                                   Vec3(x=0.0, y=0.0, z=estimated_new_side_length))
        simulation_nvt_equilibration.minimizeEnergy(tolerance=0.1 * kilojoule / mole, maxIterations=500000)
        simulation_nvt_equilibration.context.setVelocitiesToTemperature(temperature * kelvin)
        simulation_nvt_equilibration.reporters.append(StateDataReporter(f'{folder}/nvt_equilibration.csv',
                                                                        reporting_frequency_state_nvt_equilibration,
                                                                        time=True, 
                                                                        totalEnergy=True, 
                                                                        density=True, 
                                                                        volume=True,  
                                                                        potentialEnergy=True, 
                                                                        kineticEnergy=True, 
                                                                        temperature=True))
        simulation_nvt_equilibration.step(nvt_equilibration_steps)
        simulation_nvt_equilibration.saveState(f'{folder}/nvt_equilibration_state.xml')
        endpositions_nvt_equilibration = simulation_nvt_equilibration.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(topology_nvt_equilibration.topology, endpositions_nvt_equilibration, open(f'{folder}/nvt_equilibration_end.pdb', 'w'))
        df = pd.read_csv(f'{folder}/nvt_equilibration.csv')
        df.columns = ["Time (ps)","Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)","Temperature (K)","Box Volume (nm^3)","Density (g/mL)"]
        fig, axes = plt.subplots(2, 3, figsize=(40, 20))
        sns.scatterplot(data=df, ax=axes[0,0], x="Time (ps)", y="Box Volume (nm^3)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[0,1], x="Time (ps)", y="Density (g/mL)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[0,2], x="Time (ps)", y="Total Energy (kJ/mole)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[1,0], x="Time (ps)", y="Temperature (K)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[1,1], x="Time (ps)", y="Potential Energy (kJ/mole)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[1,2], x="Time (ps)", y="Kinetic Energy (kJ/mole)", s=2, marker='o')
        plt.savefig(f'{folder}/nvt_equilibration.svg', dpi=300)
        plt.savefig(f'{folder}/nvt_equilibration.png', dpi=300)

    def nve_production(self):
        """Runs a NVE production of the system during which two trajectories are written, one with unwrapped and one with wrapped coordinates.
        """

        folder = self.folder
        nvt_equilibration_state = f'{folder}/nvt_equilibration_state.xml'
        nvt_equilibration_end_simulationbox = f'{folder}/nvt_equilibration_end.pdb'
        openmm_properties = self.openmm_properties
        openmm_platform = self.openmm_platform
        pme_error_tolerance = self.pme_error_tolerance
        cutoff_distance_nm = self.cutoff_distance_nm
        cutoff_switch_distance_nm = self.cutoff_switch_distance_nm
        openmm_forcefield = self.openmm_forcefield
        production_timestep_ps = self.nve_production_timestep_fs * (10 ** (-3))
        production_steps = self.nve_production_steps
        reporting_frequency_coordinates_unwrapped = self.reporting_frequency_coordinates_unwrapped
        reporting_frequency_coordinates_wrapped = self.reporting_frequency_coordinates_wrapped
        reporting_frequency_state_nve_production = self.reporting_frequency_state_nve_production
        constraint_tolerance = self.constraint_tolerance

        topology_nve_production = PDBFile(f'{nvt_equilibration_end_simulationbox}')      
        system_nve_production = openmm_forcefield.createSystem(topology_nve_production.topology, 
                                                           nonbondedMethod=PME, 
                                                           nonbondedCutoff=cutoff_distance_nm * nanometer,
                                                           switchDistance=cutoff_switch_distance_nm * nanometer,
                                                           constraints=HBonds,
                                                           rigidWater=True,
                                                           ewaldErrorTolerance=pme_error_tolerance,
                                                           useDispersionCorrection=True)
        integrator_nve_production = VerletIntegrator(production_timestep_ps * picoseconds)
        integrator_nve_production.setConstraintTolerance(constraint_tolerance)
        simulation_nve_production = Simulation(topology_nve_production.topology, 
                                           system_nve_production, 
                                           integrator_nve_production, 
                                           openmm_platform, 
                                           openmm_properties,
                                           state=nvt_equilibration_state)
        simulation_nve_production.reporters.append(StateDataReporter(f'{folder}/nve_production.csv',
                                                                     reporting_frequency_state_nve_production,
                                                                     time=True, 
                                                                     totalEnergy=True, 
                                                                     density=True, 
                                                                     volume=True,  
                                                                     potentialEnergy=True, 
                                                                     kineticEnergy=True, 
                                                                     temperature=True))
        simulation_nve_production.reporters.append(DCDReporter(f'{folder}/nve_production_unwrapped.dcd', reporting_frequency_coordinates_unwrapped, enforcePeriodicBox=False))
        simulation_nve_production.reporters.append(DCDReporter(f'{folder}/nve_production_wrapped.dcd', reporting_frequency_coordinates_wrapped))
        simulation_nve_production.step(production_steps)
        simulation_nve_production.saveState(f'{folder}/nve_production_state.xml')
        endpositions_nve_production = simulation_nve_production.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(topology_nve_production.topology, endpositions_nve_production, open(f'{folder}/nve_production_end.pdb', 'w'))
        df = pd.read_csv(f'{folder}/nve_production.csv')
        df.columns = ["Time (ps)","Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)","Temperature (K)","Box Volume (nm^3)","Density (g/mL)"]
        fig, axes = plt.subplots(2, 3, figsize=(40, 20))
        sns.scatterplot(data=df, ax=axes[0,0], x="Time (ps)", y="Box Volume (nm^3)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[0,1], x="Time (ps)", y="Density (g/mL)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[0,2], x="Time (ps)", y="Total Energy (kJ/mole)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[1,0], x="Time (ps)", y="Temperature (K)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[1,1], x="Time (ps)", y="Potential Energy (kJ/mole)", s=2, marker='o')
        sns.scatterplot(data=df, ax=axes[1,2], x="Time (ps)", y="Kinetic Energy (kJ/mole)", s=2, marker='o')
        plt.savefig(f'{folder}/nve_production.svg', dpi=300)
        plt.savefig(f'{folder}/nve_production.png', dpi=300)