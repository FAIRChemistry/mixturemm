"""Script for simulation on bwUniCluster 2.0 and HoreKa.

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

from openmm.app import *
from openmm import *
from openmm.unit import *
import openmm.app.element as elem
import pandas as pd 
import sys
import os
from more_itertools import unique_everseen
import mdtraj as md

# get arguments from the command-line this file was started with
workspace = sys.argv[1]
simulation_properties_ = sys.argv[2]
simulation_properties = dict(x.split('=') for x in simulation_properties_.split(','))
simulation_platform = sys.argv[3]
forcefields_ = sys.argv[4]
forcefields = forcefields_.split(',')
pme_error_tolerance = float(sys.argv[5])
cutoff_distance_nm = float(sys.argv[6])
cutoff_switch_distance_nm = float(sys.argv[7])
name = sys.argv[8]
replica_number = int(sys.argv[9])
temperature = float(sys.argv[10])
npt_equilibration_subfolder = f'{workspace}/{name}/npt_equilibration'
npt_equilibration_end_simulationbox = f'{npt_equilibration_subfolder}/npt_equilibration_end.pdb'
npt_equilibration_state = f'{npt_equilibration_subfolder}/npt_equilibration_state.xml'
folder = f'{workspace}/{name}/replica_{replica_number}'
half_npt_equilibration_csv_columns = int(sys.argv[11])
nvt_equilibration_temperature_coupling_frequency = float(sys.argv[12])
nvt_equilibration_steps = int(sys.argv[13])
nvt_equilibration_timestep_ps = float(sys.argv[14]) * (10 ** (-3))
reporting_frequency_state_nvt_equilibration = int(sys.argv[15])
nvt_equilibration_state = f'{workspace}/{name}/replica_{replica_number}/nvt_equilibration_state.xml'
nvt_equilibration_end_simulationbox = f'{workspace}/{name}/replica_{replica_number}/nvt_equilibration_end.pdb'
production_timestep_ps = float(sys.argv[16]) * (10 ** (-3))
production_steps = int(sys.argv[17])
reporting_frequency_coordinates_unwrapped = int(sys.argv[18])
reporting_frequency_coordinates_wrapped = int(sys.argv[19])
reporting_frequency_state_nve_production = int(sys.argv[20])
checkpoint_frequency = int(sys.argv[21])
constraint_tolerance = float(sys.argv[22])
elements_ = sys.argv[23]
if elements_ != 'None':
    elements = dict(x.split('=') for x in elements_.split(','))
    for name, mass_ in elements.items():
        mass = float(mass_)
        _ = elem.Element(number=0, name=name, symbol=name, mass=mass*amu)
openmm_platform = Platform.getPlatformByName(f'{simulation_platform}')
openmm_forcefield = ForceField(*forcefields)

def resume_replica(operating_file, 
                   timestep, 
                   reporting_frequency_state, 
                   checkpoint_frequency, steps, 
                   traj=False, traj_pdb=None, 
                   reporting_frequency_coordinates_unwrapped=None, 
                   reporting_frequency_coordinates_wrapped=None):
    """Resumes a replica by computing and removing the csv and trajectory overhang after the last written checkpoint and returning the remaining simulation steps.

    Args:
        operating_file (str): The path to the csv file without the data ending '.csv'.
        timestep (float): The integration time step during the NVT equilibration or NVE production in femtoseconds.
        reporting_frequency_state (int): The reporting frequency of the openMM StateDataReporter in simulation steps during the NVT equilibration or NVE production.
        checkpoint_frequency (int): The frequency at which checkpoints are saved to restart from once a job continues.
        steps (int): The simulation steps of the NVT equilibration or NVE production.
        traj (bool, optional): True if trajectories are stored during the ensemble. Defaults to False.
        traj_pdb (str, optional): Path to the corresponding PDB file of the trajectory. Defaults to None.
        reporting_frequency_coordinates_unwrapped (int, optional): The reporting frequency of the trajectory with unwrapped coordinates in simulation steps during the NVE production. Defaults to None.
        reporting_frequency_coordinates_wrapped (int, optional): The reporting frequency of the trajectory with wrapped coordinates in simulation steps during the NVE production. Defaults to None.

    Returns:
        int: The remaining steps.
    """

    with open(f'{operating_file}_t.csv', 'w') as file:
        file.write('#"Time (ps)","Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)","Temperature (K)","Box Volume (nm^3)","Density (g/mL)","Elapsed Time (s)"\n')
    progress_df = pd.read_csv(f'{operating_file}.csv')
    first_line = progress_df.head(n=1)
    last_line = progress_df.tail(n=1)
    first_line.reset_index(inplace=True)
    last_line.reset_index(inplace=True)
    start_time_ = first_line.loc[0].at['#"Time (ps)"']
    progress_time_ = last_line.loc[0].at['#"Time (ps)"']
    start_time = float(start_time_)
    progress_time = float(progress_time_)
    start_steps = int(start_time / timestep)
    progress_steps = int(progress_time / timestep)
    overshoot_steps = int (progress_steps % checkpoint_frequency)
    overshoot_columns = int(overshoot_steps / reporting_frequency_state)
    progress_df.drop(progress_df.tail(n=overshoot_columns).index, inplace=True)
    progress_df.to_csv(f'{operating_file}_t.csv', index=False, header=False, mode='a')
    os.remove(f'{operating_file}.csv')
    os.rename(f'{operating_file}_t.csv', f'{operating_file}.csv')
    actual_steps = int(progress_steps - overshoot_steps)
    if traj:
        actual_frames_unwrapped = int((actual_steps - start_steps) / reporting_frequency_coordinates_unwrapped)
        tra = md.load(f'{operating_file}_unwrapped.dcd', top=f'{traj_pdb}')
        tra_new = tra[0:actual_frames_unwrapped]
        tra_new.save_dcd(f'{operating_file}_unwrapped.dcd')
        actual_frames_wrapped = int((actual_steps - start_steps) / reporting_frequency_coordinates_wrapped)
        tra = md.load(f'{operating_file}_wrapped.dcd', top=f'{traj_pdb}')
        tra_new = tra[0:actual_frames_wrapped]
        tra_new.save_dcd(f'{operating_file}_wrapped.dcd')
    new_steps = int((start_steps + steps) - (actual_steps + reporting_frequency_state))
    return new_steps

def check_if_replica_in_progress(nvt_equilibration_csv_path, nve_production_csv_path):
    """Checks if a replica simulation was already started or not yet.

    Args:
        nvt_equilibration_csv_path (str): Path to the NVT equilibration csv file.
        nve_production_csv_path (str): Path to the NVE production csv file.

    Returns:
        int: 0 if not started yet, 1 if NVT equilibration in progress and 2 if NVE production in progress.
    """

    skip = 0
    if os.path.isfile(f'{nvt_equilibration_csv_path}') == True:
        skip = 1
    if os.path.isfile(f'{nve_production_csv_path}') == True:
        skip = 2
    return skip

def nvt_equilibration_from_beginning(npt_equilibration_end_simulationbox,
                                     openmm_forcefield,
                                     cutoff_distance_nm,
                                     cutoff_switch_distance_nm,
                                     pme_error_tolerance,
                                     temperature,
                                     nvt_equilibration_temperature_coupling_frequency,
                                     nvt_equilibration_timestep_ps,
                                     constraint_tolerance,
                                     openmm_platform,
                                     simulation_properties,
                                     npt_equilibration_state,
                                     folder,
                                     npt_equilibration_subfolder,
                                     reporting_frequency_state_nvt_equilibration,
                                     checkpoint_frequency,
                                     nvt_equilibration_steps):
    """Runs a NVT equilibration from the beginning on.

    Args:
        npt_equilibration_end_simulationbox (str): Path to the simulation box PDB of the system.
        openmm_forcefield (obj): The OpenMM forcefield used to compute the forces.
        cutoff_distance_nm (float): The cutoff distance for nonbonded interactions in nanometers.
        cutoff_switch_distance_nm (_type_): Starting point of the switching function that makes the energy go smoothly to 0 at the cutoff distance in nanometers.
        pme_error_tolerance (_type_): Decides how large the grid for PME is together with the cutoff.
        temperature (float): Temperature in kelvin at which the box should be simulated.
        nvt_equilibration_temperature_coupling_frequency (float): The friction coefficient of the Langevin thermostat during the NVT equilibration in inverse picoseconds.
        nvt_equilibration_timestep_ps (float): The integration time step during the NVT equilibration in picoseconds.
        constraint_tolerance (float): The constraint tolerance specifies the distance tolerance within which constraints are maintained, as a fraction of the constrained distance.
        openmm_platform (obj): The platform the simulation should run on.
        simulation_properties (dict): The platform properties for the simulation.
        npt_equilibration_state (str): Path to the NpT equilibration end state.
        folder (str): Path to the replica folder.
        npt_equilibration_subfolder (_type_): Path to the system folder.
        reporting_frequency_state_nvt_equilibration (_type_): The reporting frequency of the openMM StateDataReporter in simulation steps during the NVT equilibration.
        checkpoint_frequency (_type_): The frequency at which checkpoints are saved to restart from once a job continues.
        nvt_equilibration_steps (_type_): The simulation steps of the NVT equilibration.
    """

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
    barostat_monte_carlo = MonteCarloBarostat(0, 0, 0)
    system_nvt_equilibration.addForce(barostat_monte_carlo)
    simulation_nvt_equilibration = Simulation(topology_nvt_equilibration.topology, 
                                              system_nvt_equilibration, 
                                              integrator_nvt_equilibration, 
                                              openmm_platform, 
                                              simulation_properties,
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
                                                                    temperature=True,
                                                                    elapsedTime=True))
    simulation_nvt_equilibration.reporters.append(CheckpointReporter(f'{folder}/nvt_equilibration.chk', checkpoint_frequency))
    simulation_nvt_equilibration.step(nvt_equilibration_steps)
    simulation_nvt_equilibration.saveState(f'{folder}/nvt_equilibration_state.xml')
    endpositions_nvt_equilibration = simulation_nvt_equilibration.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(topology_nvt_equilibration.topology, endpositions_nvt_equilibration, open(f'{folder}/nvt_equilibration_end.pdb', 'w'))

def nvt_equilibration_from_progress(npt_equilibration_end_simulationbox,
                                    openmm_forcefield,
                                    cutoff_distance_nm,
                                    cutoff_switch_distance_nm,
                                    pme_error_tolerance,
                                    temperature,
                                    nvt_equilibration_temperature_coupling_frequency,
                                    nvt_equilibration_timestep_ps,
                                    constraint_tolerance,
                                    openmm_platform,
                                    simulation_properties,
                                    folder,
                                    reporting_frequency_state_nvt_equilibration,
                                    checkpoint_frequency,
                                    nvt_equilibration_steps):
    """Runs a NVT equilibration from a progress checkpoint on.

    Args:
        npt_equilibration_end_simulationbox (str): Path to the simulation box PDB of the system.
        openmm_forcefield (obj): The OpenMM forcefield used to compute the forces.
        cutoff_distance_nm (float): The cutoff distance for nonbonded interactions in nanometers.
        cutoff_switch_distance_nm (_type_): Starting point of the switching function that makes the energy go smoothly to 0 at the cutoff distance in nanometers.
        pme_error_tolerance (_type_): Decides how large the grid for PME is together with the cutoff.
        temperature (float): Temperature in kelvin at which the box should be simulated.
        nvt_equilibration_temperature_coupling_frequency (float): The friction coefficient of the Langevin thermostat during the NVT equilibration in inverse picoseconds.
        nvt_equilibration_timestep_ps (float): The integration time step during the NVT equilibration in picoseconds.
        constraint_tolerance (float): The constraint tolerance specifies the distance tolerance within which constraints are maintained, as a fraction of the constrained distance.
        openmm_platform (obj): The platform the simulation should run on.
        simulation_properties (dict): The platform properties for the simulation.
        folder (str): Path to the replica folder.
        reporting_frequency_state_nvt_equilibration (_type_): The reporting frequency of the openMM StateDataReporter in simulation steps during the NVT equilibration.
        checkpoint_frequency (_type_): The frequency at which checkpoints are saved to restart from once a job continues.
        nvt_equilibration_steps (_type_): The simulation steps of the NVT equilibration.
    """

    new_steps = resume_replica(f'{folder}/nvt_equilibration', nvt_equilibration_timestep_ps, reporting_frequency_state_nvt_equilibration, checkpoint_frequency, nvt_equilibration_steps)
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
    barostat_monte_carlo = MonteCarloBarostat(0, 0, 0)
    system_nvt_equilibration.addForce(barostat_monte_carlo)
    simulation_nvt_equilibration = Simulation(topology_nvt_equilibration.topology, 
                                              system_nvt_equilibration, 
                                              integrator_nvt_equilibration, 
                                              openmm_platform, 
                                              simulation_properties)
    simulation_nvt_equilibration.loadCheckpoint(f'{folder}/nvt_equilibration.chk')
    datafile = open(f'{folder}/nvt_equilibration.csv', 'a')
    simulation_nvt_equilibration.reporters.append(StateDataReporter(datafile,
                                                                    reporting_frequency_state_nvt_equilibration,
                                                                    time=True, 
                                                                    totalEnergy=True, 
                                                                    density=True, 
                                                                    volume=True,  
                                                                    potentialEnergy=True, 
                                                                    kineticEnergy=True, 
                                                                    temperature=True,
                                                                    elapsedTime=True))
    simulation_nvt_equilibration.reporters.append(CheckpointReporter(f'{folder}/nvt_equilibration.chk', checkpoint_frequency))
    simulation_nvt_equilibration.step(new_steps)
    simulation_nvt_equilibration.saveState(f'{folder}/nvt_equilibration_state.xml')
    endpositions_nvt_equilibration = simulation_nvt_equilibration.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(topology_nvt_equilibration.topology, endpositions_nvt_equilibration, open(f'{folder}/nvt_equilibration_end.pdb', 'w'))
    datafile.close()
    with open(f'{folder}/nvt_equilibration.csv','r') as f, open(f'{folder}/nvt_equilibration_t.csv','w') as out_file:
        out_file.writelines(unique_everseen(f))
    os.remove(f'{folder}/nvt_equilibration.csv')
    os.rename(f'{folder}/nvt_equilibration_t.csv', f'{folder}/nvt_equilibration.csv')

def nve_production_from_beginning(nvt_equilibration_end_simulationbox,
                                  openmm_forcefield,
                                  cutoff_distance_nm,
                                  cutoff_switch_distance_nm,
                                  pme_error_tolerance,
                                  production_timestep_ps,
                                  constraint_tolerance,
                                  openmm_platform,
                                  simulation_properties,
                                  nvt_equilibration_state,
                                  folder,
                                  reporting_frequency_state_nve_production,
                                  reporting_frequency_coordinates_unwrapped,
                                  reporting_frequency_coordinates_wrapped,
                                  checkpoint_frequency,
                                  production_steps):
    """Runs a NVE production from the beginning on.

    Args:
        nvt_equilibration_end_simulationbox (str): Path to the simulation box PDB of the end state of the NVT equilibration.
        openmm_forcefield (obj): The OpenMM forcefield used to compute the forces.
        cutoff_distance_nm (float): The cutoff distance for nonbonded interactions in nanometers.
        cutoff_switch_distance_nm (float): Starting point of the switching function that makes the energy go smoothly to 0 at the cutoff distance in nanometers.
        pme_error_tolerance (float): Decides how large the grid for PME is together with the cutoff.
        production_timestep_ps (float): The integration time step during the NVE production in picoseconds.
        constraint_tolerance (float): The constraint tolerance specifies the distance tolerance within which constraints are maintained, as a fraction of the constrained distance.
        openmm_platform (obj): The platform the simulation should run on.
        simulation_properties (dict): The platform properties for the simulation.
        nvt_equilibration_state (str): Path to the NVT equilibration end state.
        folder (str): Path to the replica folder.
        reporting_frequency_state_nve_production (int): The reporting frequency of the openMM StateDataReporter in simulation steps during the NVE production.
        reporting_frequency_coordinates_unwrapped (int): The reporting frequency of the trajectory with unwrapped coordinates in simulation steps during the NVE production.
        reporting_frequency_coordinates_wrapped (int): The reporting frequency of the trajectory with wrapped coordinates in simulation steps during the NVE production.
        checkpoint_frequency (int): The frequency at which checkpoints are saved to restart from once a job continues.
        production_steps (int): The simulation steps of the NVE production.
    """

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
    barostat_monte_carlo = MonteCarloBarostat(0, 0, 0)
    system_nve_production.addForce(barostat_monte_carlo)
    simulation_nve_production = Simulation(topology_nve_production.topology, 
                                           system_nve_production, 
                                           integrator_nve_production, 
                                           openmm_platform, 
                                           simulation_properties,
                                           state=nvt_equilibration_state)
    simulation_nve_production.reporters.append(StateDataReporter(f'{folder}/nve_production.csv',
                                                                reporting_frequency_state_nve_production,
                                                                time=True, 
                                                                totalEnergy=True, 
                                                                density=True, 
                                                                volume=True,  
                                                                potentialEnergy=True, 
                                                                kineticEnergy=True, 
                                                                temperature=True,
                                                                elapsedTime=True))
    simulation_nve_production.reporters.append(DCDReporter(f'{folder}/nve_production_unwrapped.dcd', reporting_frequency_coordinates_unwrapped, enforcePeriodicBox=False))
    simulation_nve_production.reporters.append(DCDReporter(f'{folder}/nve_production_wrapped.dcd', reporting_frequency_coordinates_wrapped))
    simulation_nve_production.reporters.append(CheckpointReporter(f'{folder}/nve_production.chk', checkpoint_frequency))
    simulation_nve_production.step(production_steps)
    simulation_nve_production.saveState(f'{folder}/nve_production_state.xml')
    endpositions_nve_production = simulation_nve_production.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(topology_nve_production.topology, endpositions_nve_production, open(f'{folder}/nve_production_end.pdb', 'w'))

def nve_production_from_progress(nvt_equilibration_end_simulationbox,
                                 openmm_forcefield,
                                 cutoff_distance_nm,
                                 cutoff_switch_distance_nm,
                                 pme_error_tolerance,
                                 production_timestep_ps,
                                 constraint_tolerance,
                                 openmm_platform,
                                 simulation_properties,
                                 folder,
                                 reporting_frequency_state_nve_production,
                                 reporting_frequency_coordinates_unwrapped,
                                 reporting_frequency_coordinates_wrapped,
                                 checkpoint_frequency,
                                 production_steps):
    """Runs a NVE production from a progress checkpoint on.

    Args:
        nvt_equilibration_end_simulationbox (str): Path to the simulation box PDB of the end state of the NVT equilibration.
        openmm_forcefield (obj): The OpenMM forcefield used to compute the forces.
        cutoff_distance_nm (float): The cutoff distance for nonbonded interactions in nanometers.
        cutoff_switch_distance_nm (float): Starting point of the switching function that makes the energy go smoothly to 0 at the cutoff distance in nanometers.
        pme_error_tolerance (float): Decides how large the grid for PME is together with the cutoff.
        production_timestep_ps (float): The integration time step during the NVE production in picoseconds.
        constraint_tolerance (float): The constraint tolerance specifies the distance tolerance within which constraints are maintained, as a fraction of the constrained distance.
        openmm_platform (obj): The platform the simulation should run on.
        simulation_properties (dict): The platform properties for the simulation.
        folder (str): Path to the replica folder.
        reporting_frequency_state_nve_production (int): The reporting frequency of the openMM StateDataReporter in simulation steps during the NVE production.
        reporting_frequency_coordinates_unwrapped (int): The reporting frequency of the trajectory with unwrapped coordinates in simulation steps during the NVE production.
        reporting_frequency_coordinates_wrapped (int): The reporting frequency of the trajectory with wrapped coordinates in simulation steps during the NVE production.
        checkpoint_frequency (int): The frequency at which checkpoints are saved to restart from once a job continues.
        production_steps (int): The simulation steps of the NVE production.
    """

    new_steps = resume_replica(f'{folder}/nve_production', production_timestep_ps, reporting_frequency_state_nve_production, checkpoint_frequency, production_steps, traj=True, traj_pdb=nvt_equilibration_end_simulationbox, reporting_frequency_coordinates_unwrapped=reporting_frequency_coordinates_unwrapped, reporting_frequency_coordinates_wrapped=reporting_frequency_coordinates_wrapped)
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
    barostat_monte_carlo = MonteCarloBarostat(0, 0, 0)
    system_nve_production.addForce(barostat_monte_carlo)
    simulation_nve_production = Simulation(topology_nve_production.topology, 
                                           system_nve_production, 
                                           integrator_nve_production, 
                                           openmm_platform, 
                                           simulation_properties)
    simulation_nve_production.loadCheckpoint(f'{folder}/nve_production.chk')
    datafile = open(f'{folder}/nve_production.csv', 'a')
    simulation_nve_production.reporters.append(StateDataReporter(datafile,
                                                                reporting_frequency_state_nve_production,
                                                                time=True, 
                                                                totalEnergy=True, 
                                                                density=True, 
                                                                volume=True,  
                                                                potentialEnergy=True, 
                                                                kineticEnergy=True, 
                                                                temperature=True,
                                                                elapsedTime=True))
    simulation_nve_production.reporters.append(DCDReporter(f'{folder}/nve_production_unwrapped.dcd', reporting_frequency_coordinates_unwrapped, append=True, enforcePeriodicBox=False))
    simulation_nve_production.reporters.append(DCDReporter(f'{folder}/nve_production_wrapped.dcd', reporting_frequency_coordinates_wrapped, append=True))
    simulation_nve_production.reporters.append(CheckpointReporter(f'{folder}/nve_production.chk', checkpoint_frequency))
    simulation_nve_production.step(new_steps)
    simulation_nve_production.saveState(f'{folder}/nve_production_state.xml')
    endpositions_nve_production = simulation_nve_production.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(topology_nve_production.topology, endpositions_nve_production, open(f'{folder}/nve_production_end.pdb', 'w'))
    datafile.close()
    with open(f'{folder}/nve_production.csv','r') as f, open(f'{folder}/nve_production_t.csv','w') as out_file:
        out_file.writelines(unique_everseen(f))
    os.remove(f'{folder}/nve_production.csv')
    os.rename(f'{folder}/nve_production_t.csv', f'{folder}/nve_production.csv')

# This is the part that runs
skip = check_if_replica_in_progress(f'{folder}/nvt_equilibration.csv', f'{folder}/nve_production.csv')
if skip == 0:
    nvt_equilibration_from_beginning(npt_equilibration_end_simulationbox,
                                     openmm_forcefield,
                                     cutoff_distance_nm,
                                     cutoff_switch_distance_nm,
                                     pme_error_tolerance,
                                     temperature,
                                     nvt_equilibration_temperature_coupling_frequency,
                                     nvt_equilibration_timestep_ps,
                                     constraint_tolerance,
                                     openmm_platform,
                                     simulation_properties,
                                     npt_equilibration_state,
                                     folder,
                                     npt_equilibration_subfolder,
                                     reporting_frequency_state_nvt_equilibration,
                                     checkpoint_frequency,
                                     nvt_equilibration_steps)
if skip == 1:
    nvt_equilibration_from_progress(npt_equilibration_end_simulationbox,
                                    openmm_forcefield,
                                    cutoff_distance_nm,
                                    cutoff_switch_distance_nm,
                                    pme_error_tolerance,
                                    temperature,
                                    nvt_equilibration_temperature_coupling_frequency,
                                    nvt_equilibration_timestep_ps,
                                    constraint_tolerance,
                                    openmm_platform,
                                    simulation_properties,
                                    folder,
                                    reporting_frequency_state_nvt_equilibration,
                                    checkpoint_frequency,
                                    nvt_equilibration_steps)
if skip < 2:
    nve_production_from_beginning(nvt_equilibration_end_simulationbox,
                                  openmm_forcefield,
                                  cutoff_distance_nm,
                                  cutoff_switch_distance_nm,
                                  pme_error_tolerance,
                                  production_timestep_ps,
                                  constraint_tolerance,
                                  openmm_platform,
                                  simulation_properties,
                                  nvt_equilibration_state,
                                  folder,
                                  reporting_frequency_state_nve_production,
                                  reporting_frequency_coordinates_unwrapped,
                                  reporting_frequency_coordinates_wrapped,
                                  checkpoint_frequency,
                                  production_steps)
if skip == 2:
    nve_production_from_progress(nvt_equilibration_end_simulationbox,
                                 openmm_forcefield,
                                 cutoff_distance_nm,
                                 cutoff_switch_distance_nm,
                                 pme_error_tolerance,
                                 production_timestep_ps,
                                 constraint_tolerance,
                                 openmm_platform,
                                 simulation_properties,
                                 folder,
                                 reporting_frequency_state_nve_production,
                                 reporting_frequency_coordinates_unwrapped,
                                 reporting_frequency_coordinates_wrapped,
                                 checkpoint_frequency,
                                 production_steps)

import seaborn as sns
import matplotlib.pyplot as plt
plt.switch_backend('agg')
plt.ioff()

df = pd.read_csv(f'{folder}/nvt_equilibration.csv')
df.columns = ["Time (ps)","Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)","Temperature (K)","Box Volume (nm^3)","Density (g/mL)","Elapsed Time (s)"]
fig, axes = plt.subplots(2, 3, figsize=(40, 20))
sns.scatterplot(data=df, ax=axes[0,0], x="Time (ps)", y="Box Volume (nm^3)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[0,1], x="Time (ps)", y="Density (g/mL)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[0,2], x="Time (ps)", y="Total Energy (kJ/mole)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[1,0], x="Time (ps)", y="Temperature (K)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[1,1], x="Time (ps)", y="Potential Energy (kJ/mole)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[1,2], x="Time (ps)", y="Kinetic Energy (kJ/mole)", s=2, marker='o')
plt.savefig(f'{folder}/nvt_equilibration.svg', dpi=300)
plt.savefig(f'{folder}/nvt_equilibration.png', dpi=300)
plt.clf()

df = pd.read_csv(f'{folder}/nve_production.csv')
df.columns = ["Time (ps)","Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)","Temperature (K)","Box Volume (nm^3)","Density (g/mL)","Elapsed Time (s)"]
fig, axes = plt.subplots(2, 3, figsize=(40, 20))
sns.scatterplot(data=df, ax=axes[0,0], x="Time (ps)", y="Box Volume (nm^3)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[0,1], x="Time (ps)", y="Density (g/mL)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[0,2], x="Time (ps)", y="Total Energy (kJ/mole)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[1,0], x="Time (ps)", y="Temperature (K)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[1,1], x="Time (ps)", y="Potential Energy (kJ/mole)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[1,2], x="Time (ps)", y="Kinetic Energy (kJ/mole)", s=2, marker='o')
plt.savefig(f'{folder}/nve_production.svg', dpi=300)
plt.savefig(f'{folder}/nve_production.png', dpi=300)
plt.clf()