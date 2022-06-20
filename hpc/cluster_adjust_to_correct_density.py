"""Script for adjusting to correct density on bwUniCluster 2.0 and HoreKa.

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
from more_itertools import unique_everseen
import pandas as pd 
import numpy as np
import sys
import os

# get arguments from the command-line this file was started with
START_TEMPERATURE = 1
TEMPERATURE_STEP = 0.1
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
temperature = float(sys.argv[9])
pressure = float(sys.argv[10])
simulationbox_name = sys.argv[11]
folder = f'{workspace}/{name}/npt_equilibration'
simulationbox_path = f'{workspace}/boxes/mixture_{simulationbox_name}.pdb'
npt_equilibration_pressure_coupling_frequency = int(sys.argv[12])
npt_equilibration_temperature_coupling_frequency = float(sys.argv[13])
npt_equilibration_timestep_ps = float(sys.argv[14]) * (10 ** (-3))
npt_equilibration_steps = int(sys.argv[15])
reporting_frequency_state_npt_equilibration = int(sys.argv[16])
constraint_tolerance = float(sys.argv[17])
checkpoint_frequency = int(sys.argv[18])
elements_ = sys.argv[19]
if elements_ != 'None':
    elements = dict(x.split('=') for x in elements_.split(','))
    for name, mass_ in elements.items():
        mass = float(mass_)
        _ = elem.Element(number=0, name=name, symbol=name, mass=mass*amu)
heating_interval_steps = int(1000 / (npt_equilibration_timestep_ps * 10 ** 3))
openmm_platform = Platform.getPlatformByName(f'{simulation_platform}')
openmm_forcefield = ForceField(*forcefields)
heating_steps = (len(np.arange(START_TEMPERATURE, temperature, TEMPERATURE_STEP)) * heating_interval_steps + 10 * heating_interval_steps)

def resume_system(operating_file, 
                  timestep, 
                  reporting_frequency_state, 
                  checkpoint_frequency, 
                  steps, 
                  heating_steps):
    """Resumes a NpT equilibration by computing and removing the csv overhang after the last written checkpoint and returning the remaining simulation steps.

    Args:
        operating_file (str): The path to the csv file without the data ending '.csv'.
        timestep (float): The integration time step during the NpT equilibration in femtoseconds.
        reporting_frequency_state (int): The reporting frequency of the openMM StateDataReporter in simulation steps during the NpT equilibration.
        checkpoint_frequency (int): The frequency at which checkpoints are saved to restart from once a job continues.
        steps (int): The simulation steps of the NpT equilibration.
        heating_steps (int): The simulation steps to advance during one heating step.

    Returns:
        int: The remaining steps, wheter the simulation is still in the heating phase and the remaining heating steps.
    """

    with open(f'{operating_file}_t.csv', 'w') as file:
        file.write('#"Time (ps)","Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)","Temperature (K)","Box Volume (nm^3)","Density (g/mL)","Elapsed Time (s)"\n')
    progress_df = pd.read_csv(f'{operating_file}.csv')
    first_line = progress_df.head(n=1)
    last_line = progress_df.tail(n=1)
    first_line.reset_index(inplace=True)
    last_line.reset_index(inplace=True)
    start_time = first_line.loc[0].at['#"Time (ps)"']
    progress_time = last_line.loc[0].at['#"Time (ps)"']
    start_steps = int(start_time / timestep)
    progress_steps = int(progress_time / timestep)
    overshoot_steps = int (progress_steps % checkpoint_frequency)
    overshoot_columns = int(overshoot_steps / reporting_frequency_state)
    progress_df.drop(progress_df.tail(n=overshoot_columns).index, inplace=True)
    progress_df.to_csv(f'{operating_file}_t.csv', index=False, header=False, mode='a')
    os.remove(f'{operating_file}.csv')
    os.rename(f'{operating_file}_t.csv', f'{operating_file}.csv')
    actual_steps = int(progress_steps - overshoot_steps)
    heating = False
    if actual_steps < heating_steps:
        heating = True
        new_steps = steps
        heating_portion = int((start_steps + heating_steps) - (actual_steps + reporting_frequency_state))
    else:
        new_steps = int((start_steps + steps + heating_steps) - (actual_steps + reporting_frequency_state))
    return new_steps, heating, heating_portion

def check_if_system_in_progress(npt_equilibration_csv_path):
    """Checks if a system simulation was already started or not yet.

    Args:
        npt_equilibration_csv_path (str): Path to the NpT equilibration csv file.

    Returns:
        int: 0 if not started yet and 1 if in progress.
    """

    skip = 0
    if os.path.isfile(f'{npt_equilibration_csv_path}') == True:
        skip = 1
    return skip

def npt_equilibration_from_beginning(simulationbox_path, 
                                     openmm_forcefield,
                                     cutoff_distance_nm,
                                     cutoff_switch_distance_nm,
                                     pme_error_tolerance, 
                                     START_TEMPERATURE, 
                                     npt_equilibration_timestep_ps,
                                     constraint_tolerance,
                                     pressure,
                                     temperature,
                                     openmm_platform,
                                     simulation_properties,
                                     folder,
                                     reporting_frequency_state_npt_equilibration,
                                     checkpoint_frequency,
                                     TEMPERATURE_STEP,
                                     heating_interval_steps,
                                     npt_equilibration_temperature_coupling_frequency,
                                     npt_equilibration_pressure_coupling_frequency,
                                     npt_equilibration_steps):
    """Runs a NpT equilibration from the beginning on.

    Args:
        simulationbox_path (str): Path to the simulation box PDB of the system.
        openmm_forcefield (obj): The OpenMM forcefield used to compute the forces.
        cutoff_distance_nm (float): The cutoff distance for nonbonded interactions in nanometers.
        cutoff_switch_distance_nm (float): Starting point of the switching function that makes the energy go smoothly to 0 at the cutoff distance in nanometers.
        pme_error_tolerance (float): Decides how large the grid for PME is together with the cutoff.
        START_TEMPERATURE (int): The starting temperature of the heating (1 Kelvin).
        npt_equilibration_timestep_ps (float): The integration time step during the NpT equilibration in picoseconds.
        constraint_tolerance (float): The constraint tolerance specifies the distance tolerance within which constraints are maintained, as a fraction of the constrained distance.
        pressure (float): Pressure in bar at which the box should be simulated.
        temperature (float): Temperature in kelvin at which the box should be simulated.
        openmm_platform (obj): The platform the simulation should run on.
        simulation_properties (dict): The platform properties for the simulation.
        folder (str): Path to the system folder.
        reporting_frequency_state_npt_equilibration (int): The reporting frequency of the openMM StateDataReporter in simulation steps during the NpT equilibration.
        checkpoint_frequency (int): The frequency at which checkpoints are saved to restart from once a job continues.
        TEMPERATURE_STEP (float): The temperature jump of one heating step (0.1 Kelvin).
        heating_interval_steps (int): The simulation steps to advance during one heating step.
        npt_equilibration_temperature_coupling_frequency (float): The friction coefficient of the Langevin thermostat during the NpT equilibration in inverse picoseconds.
        npt_equilibration_pressure_coupling_frequency (int): The pressure coupling frequency at which the Monte Carlo barostat should interact with the system and attempt a Monte Carlo move to adjust the volume during the NpT equilibration in simulation steps.
        npt_equilibration_steps (int): The simulation steps of the NpT equilibration.
    """

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
                                            simulation_properties)
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
                                                                    temperature=True,
                                                                    elapsedTime=True))
    simulation_npt_equilibration.reporters.append(CheckpointReporter(f'{folder}/npt_equilibration.chk', checkpoint_frequency))
    for temp in np.arange(START_TEMPERATURE, temperature, TEMPERATURE_STEP):
        integrator_npt_equilibration.setTemperature(temp * kelvin)
        simulation_npt_equilibration.step(heating_interval_steps)
    simulation_npt_equilibration.step(heating_interval_steps * 10)
    integrator_npt_equilibration.setTemperature(temperature * kelvin)
    integrator_npt_equilibration.setFriction(npt_equilibration_temperature_coupling_frequency)
    barostat_monte_carlo.setFrequency(npt_equilibration_pressure_coupling_frequency)
    simulation_npt_equilibration.step(npt_equilibration_steps)
    simulation_npt_equilibration.saveState(f'{folder}/npt_equilibration_state.xml')
    endpositions_npt_equilibration = simulation_npt_equilibration.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(topology_npt_equilibration.topology, endpositions_npt_equilibration, open(f'{folder}/npt_equilibration_end.pdb', 'w'))

def npt_equilibration_from_progress(simulationbox_path,  
                                    openmm_forcefield,
                                    pme_error_tolerance, 
                                    START_TEMPERATURE, 
                                    npt_equilibration_timestep_ps,
                                    constraint_tolerance,
                                    pressure,
                                    temperature,
                                    openmm_platform,
                                    simulation_properties,
                                    folder,
                                    reporting_frequency_state_npt_equilibration,
                                    checkpoint_frequency,
                                    TEMPERATURE_STEP,
                                    heating_interval_steps,
                                    npt_equilibration_temperature_coupling_frequency,
                                    npt_equilibration_pressure_coupling_frequency,
                                    npt_equilibration_steps,
                                    heating_steps):
    """Runs a NpT equilibration from a progress checkpoint on.

    Args:
        simulationbox_path (str): Path to the simulation box PDB of the system.
        openmm_forcefield (obj): The OpenMM forcefield used to compute the forces.
        cutoff_distance_nm (float): The cutoff distance for nonbonded interactions in nanometers.
        cutoff_switch_distance_nm (float): Starting point of the switching function that makes the energy go smoothly to 0 at the cutoff distance in nanometers.
        pme_error_tolerance (float): Decides how large the grid for PME is together with the cutoff.
        START_TEMPERATURE (int): The starting temperature of the heating (1 Kelvin).
        npt_equilibration_timestep_ps (float): The integration time step during the NpT equilibration in picoseconds.
        constraint_tolerance (float): The constraint tolerance specifies the distance tolerance within which constraints are maintained, as a fraction of the constrained distance.
        pressure (float): Pressure in bar at which the box should be simulated.
        temperature (float): Temperature in kelvin at which the box should be simulated.
        openmm_platform (obj): The platform the simulation should run on.
        simulation_properties (dict): The platform properties for the simulation.
        folder (str): Path to the system folder.
        reporting_frequency_state_npt_equilibration (int): The reporting frequency of the openMM StateDataReporter in simulation steps during the NpT equilibration.
        checkpoint_frequency (int): The frequency at which checkpoints are saved to restart from once a job continues.
        TEMPERATURE_STEP (float): The temperature jump of one heating step (0.1 Kelvin).
        heating_interval_steps (int): The simulation steps to advance during one heating step.
        npt_equilibration_temperature_coupling_frequency (float): The friction coefficient of the Langevin thermostat during the NpT equilibration in inverse picoseconds.
        npt_equilibration_pressure_coupling_frequency (int): The pressure coupling frequency at which the Monte Carlo barostat should interact with the system and attempt a Monte Carlo move to adjust the volume during the NpT equilibration in simulation steps.
        npt_equilibration_steps (int): The simulation steps of the NpT equilibration.
    """

    new_steps, heating, heating_portion = resume_system(f'{folder}/npt_equilibration', npt_equilibration_timestep_ps, reporting_frequency_state_npt_equilibration, checkpoint_frequency, npt_equilibration_steps, heating_steps)
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
                                            simulation_properties)
    simulation_npt_equilibration.loadCheckpoint(f'{folder}/npt_equilibration.chk')
    datafile = open(f'{folder}/npt_equilibration.csv', 'a')
    simulation_npt_equilibration.reporters.append(StateDataReporter(datafile,
                                                                    reporting_frequency_state_npt_equilibration,
                                                                    time=True, 
                                                                    totalEnergy=True, 
                                                                    density=True, 
                                                                    volume=True,  
                                                                    potentialEnergy=True, 
                                                                    kineticEnergy=True, 
                                                                    temperature=True,
                                                                    elapsedTime=True))
    if heating:
        temp_arr = np.arange(START_TEMPERATURE, temperature, TEMPERATURE_STEP)
        progress_ratio = heating_portion / heating_steps
        progress_slice = int(len(temp_arr) * progress_ratio)
        temp_arr_new = temp_arr[progress_slice:]
        for temp in temp_arr_new:
            integrator_npt_equilibration.setTemperature(temp * kelvin)
            simulation_npt_equilibration.step(heating_interval_steps)
        simulation_npt_equilibration.step(heating_interval_steps * 10)
    integrator_npt_equilibration.setTemperature(temperature * kelvin)
    integrator_npt_equilibration.setFriction(npt_equilibration_temperature_coupling_frequency)
    barostat_monte_carlo.setFrequency(npt_equilibration_pressure_coupling_frequency)
    simulation_npt_equilibration.step(new_steps)
    simulation_npt_equilibration.saveState(f'{folder}/npt_equilibration_state.xml')
    endpositions_npt_equilibration = simulation_npt_equilibration.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(topology_npt_equilibration.topology, endpositions_npt_equilibration, open(f'{folder}/npt_equilibration_end.pdb', 'w'))
    datafile.close()
    with open(f'{folder}/npt_equilibration.csv','r') as f, open(f'{folder}/npt_equilibration_t.csv','w') as out_file:
        out_file.writelines(unique_everseen(f))
    os.remove(f'{folder}/npt_equilibration.csv')
    os.rename(f'{folder}/npt_equilibration_t.csv', f'{folder}/npt_equilibration.csv')

# This is the part that runs
skip = check_if_system_in_progress(f'{folder}/npt_equilibration.csv')
if skip == 0:
    npt_equilibration_from_beginning(simulationbox_path, 
                                     openmm_forcefield,
                                     cutoff_distance_nm,
                                     cutoff_switch_distance_nm,
                                     pme_error_tolerance, 
                                     START_TEMPERATURE, 
                                     npt_equilibration_timestep_ps,
                                     constraint_tolerance,
                                     pressure,
                                     temperature,
                                     openmm_platform,
                                     simulation_properties,
                                     folder,
                                     reporting_frequency_state_npt_equilibration,
                                     checkpoint_frequency,
                                     TEMPERATURE_STEP,
                                     heating_interval_steps,
                                     npt_equilibration_temperature_coupling_frequency,
                                     npt_equilibration_pressure_coupling_frequency,
                                     npt_equilibration_steps)
if skip == 1:
    npt_equilibration_from_progress(simulationbox_path,  
                                    openmm_forcefield,
                                    pme_error_tolerance, 
                                    START_TEMPERATURE, 
                                    npt_equilibration_timestep_ps,
                                    constraint_tolerance,
                                    pressure,
                                    temperature,
                                    openmm_platform,
                                    simulation_properties,
                                    folder,
                                    reporting_frequency_state_npt_equilibration,
                                    checkpoint_frequency,
                                    TEMPERATURE_STEP,
                                    heating_interval_steps,
                                    npt_equilibration_temperature_coupling_frequency,
                                    npt_equilibration_pressure_coupling_frequency,
                                    npt_equilibration_steps,
                                    heating_steps)


import seaborn as sns
import matplotlib.pyplot as plt
plt.switch_backend('agg')
plt.ioff()

df = pd.read_csv(f'{folder}/npt_equilibration.csv')
df.columns = ["Time (ps)","Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)","Temperature (K)","Box Volume (nm^3)","Density (g/mL)","Elapsed Time (s)"]
fig, axes = plt.subplots(2, 3, figsize=(40, 20))
sns.scatterplot(data=df, ax=axes[0,0], x="Time (ps)", y="Box Volume (nm^3)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[0,1], x="Time (ps)", y="Density (g/mL)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[0,2], x="Time (ps)", y="Total Energy (kJ/mole)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[1,0], x="Time (ps)", y="Temperature (K)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[1,1], x="Time (ps)", y="Potential Energy (kJ/mole)", s=2, marker='o')
sns.scatterplot(data=df, ax=axes[1,2], x="Time (ps)", y="Kinetic Energy (kJ/mole)", s=2, marker='o')
plt.savefig(f'{folder}/npt_equilibration.svg', dpi=300)
plt.savefig(f'{folder}/npt_equilibration.png', dpi=300)