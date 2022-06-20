"""Script for analysis on bwUniCluster 2.0 and HoreKa.

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

import sys
import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
import csv
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
import json
plt.switch_backend('agg')
plt.ioff()

# get arguments from the command-line this file was started with
SELF_DIFFUSION_COEFFICIENT_CONVERSION_FACTOR = 1 * ((10 ** (-20)) / (10 ** (-12)))
workspace = sys.argv[1]
name = sys.argv[2]
replica_number = int(sys.argv[3])
folder = f'{workspace}/{name}/replica_{replica_number}'
resdir = f'{workspace}/results'
res_rawdir = f'{workspace}/results/raw_data'
partsdir = f'{workspace}/results/parts'
chi_water = float(sys.argv[4])
water_abbreviation = sys.argv[5]
molecule_abbreviations_ = sys.argv[6]
molecule_abbreviations = molecule_abbreviations_.split(',')
time_between_frames = float(sys.argv[7])
fit_starting_frame = int(sys.argv[8])
fit_ending_frame = int(sys.argv[9])
temperature = float(sys.argv[10])
pressure = float(sys.argv[11])
total_number_molecules = int(sys.argv[12])
molecule_names_ = sys.argv[13]
molecule_names = molecule_names_.split(',')
water_name = sys.argv[14]

def load_mdanalysis_universe(path_to_pdb, path_to_trajectory):
    """Loads a MDAnalysis universe from a PDB file and a corresponding trajectory.

    Args:
        path_to_pdb (str): Path to PDB file.
        path_to_trajectory (str): Path to trajectory.

    Returns:
        obj: MDAnalysis universe.
    """

    try:
        mdanalysis_universe = mda.Universe(f'{path_to_pdb}')
        mdanalysis_universe.load_new(f'{path_to_trajectory}')
    except:
        mdanalysis_universe = None
    return mdanalysis_universe

def compute_msd_and_write_to_csv(mdanalysis_universe, abbreviation, time_between_frames, res_rawdir, name, replica_number):
    """Computes the mean squared displacement using MDAnalysis and writes the results to a CSV file.

    Args:
        mdanalysis_universe (obj): MDAnalysis universe.
        abbreviation (str): The abbreviation of a molecule that is used inside the PDB file.
        time_between_frames (float): The time between frames in the trajectory.
        res_rawdir (str): Directory with raw, unprocessed results.
        name (str): System name.
        replica_number (int): Replica number.

    Returns:
        list: MSD time axis and corresponding values.
    """

    Msd = msd.EinsteinMSD(mdanalysis_universe, select=f'resname {abbreviation}', msd_type='xyz', fft=True)
    Msd.run()
    msd_timeseries = Msd.results.timeseries
    msd_number_frames = Msd.n_frames
    msd_time_axis = list(np.arange(msd_number_frames) * time_between_frames)
    rows_for_csv = zip(msd_time_axis, msd_timeseries)
    with open(f'{res_rawdir}/{name}_{abbreviation}_replica_{replica_number}_msd.csv', 'w') as f:
        writer = csv.writer(f)
        for row in rows_for_csv:
            writer.writerow(row)
    del rows_for_csv
    del Msd
    return msd_time_axis, msd_timeseries

def fit_objective(x, a, b):
    """Fitting objective for computing self-diffusion coefficients from mean squared displacement.

    Args:
        x (float): X-axis values.
        a (float): Slope of the linear fit.
        b (float): Y-axis intercept of the linear fit.

    Returns:
        tuple: The y-axis value of the linear fit.
    """

    return (a * x + b)

def linear_fit_self_diffusion_coefficients_and_plot_msd(msd_time_axis, msd_timeseries, res_rawdir, name, abbreviation):
    """Performs the linear fit of the mean squared displacement and plots it.

    Args:
        msd_time_axis (list): MSD time axis.
        msd_timeseries (list): Corresponding mean squared displacement of the MSD time axis.
        res_rawdir (str): Directory with raw, unprocessed results.
        name (str): System name.
        abbreviation (str): The abbreviation of a molecule that is used inside the PDB file.

    Returns:
        float: The slope of the linear fit.
    """

    try:
        popt, pcov = curve_fit(fit_objective, msd_time_axis[fit_starting_frame:fit_ending_frame], msd_timeseries[fit_starting_frame:fit_ending_frame])
        a, b = popt
        fit_slope = a
        sns.scatterplot(x=msd_time_axis[fit_starting_frame:fit_ending_frame], y=msd_timeseries[fit_starting_frame:fit_ending_frame], s=2, marker='o')
        x_line = np.arange(msd_time_axis[fit_starting_frame], msd_time_axis[fit_ending_frame], 1)
        y_line = fit_objective(x_line, a, b)
        plt.plot(x_line, y_line, '--', color='red', label='fit: a=%5.8f, b=%5.3f' %tuple(popt))
        plt.legend(loc='upper left')
        plt.savefig(f'{res_rawdir}/{name}_{abbreviation}_replica_{replica_number}_msd.png', dpi=300)
        plt.clf()
    except:
        fit_slope = None
        sns.scatterplot(x=msd_time_axis[fit_starting_frame:fit_ending_frame], y=msd_timeseries[fit_starting_frame:fit_ending_frame], s=2, marker='o')
        plt.savefig(f'{res_rawdir}/{name}_{abbreviation}_replica_{replica_number}_msd.png', dpi=300)
        plt.clf()
    return fit_slope

def write_self_diffusion_coefficients_to_json(resdir, 
                                              temperature, 
                                              pressure, 
                                              chi_water, 
                                              total_number_molecules, 
                                              fit_slope, 
                                              folder, 
                                              replica_number, 
                                              molecule_names,
                                              molecule_abbreviations,
                                              abbreviation,
                                              partsdir,
                                              name):
    """Writes self-diffusion coefficients part JSON files containing the self-diffusion coefficient one molecule of a replica.

    Args:
        resdir (str): Path to the results directory.
        temperature (float): Temperature in kelvin at which the box was simulated.
        pressure (float): Pressure in bar at which the box was simulated.
        chi_water (float): Water mole fraction of the box.
        total_number_molecules (int): Total number of molecules inside the simulation box.
        fit_slope (float): The slope of the linear MSD fit.
        folder (str): Path to replica folder.
        replica_number (int): Replica number.
        molecule_names (list):  A list with all molecule names of the mixture.
        molecule_abbreviations (list): A list with all molecule abbreviations of the mixture.
        abbreviation (str): The abbreviation of a molecule that is used inside the PDB file.
        partsdir (str): Path to the parts directory.
        name (str): System name.
    """

    with open(f'{resdir}/densities.json', 'r', encoding='utf-8') as f:
        data = json.load(f)
        combination = (temperature, pressure, chi_water, total_number_molecules)
        for entry in data['densities']:
            if combination == (entry['temperature'], entry['pressure'], entry['chi_water'], entry['total_number_molecules']):
                volume = entry['volume']              
        if fit_slope is None:
            with open(f'{resdir}/msd_fit_error_report.txt', 'a+') as f:
                f.write(f'{folder} could not be fitted.')
        else:
            self_diffusion_coefficient_raw = fit_slope * (1 / 6)
            self_diffusion_coefficient = self_diffusion_coefficient_raw * SELF_DIFFUSION_COEFFICIENT_CONVERSION_FACTOR
            json_entry = {
                'total_number_molecules' : total_number_molecules,
                'replica_number' : replica_number,
                'temperature' : temperature,
                'pressure' : pressure,
                'chi_water' : chi_water,
                'volume' : volume,
                'molecules' : molecule_names,
                'abbreviations' : molecule_abbreviations,
                'abbreviation' : abbreviation,
                'self_diffusion_coefficient' : self_diffusion_coefficient
            }
            with open(f'{partsdir}/{name}_{abbreviation}_replica_{replica_number}_sdc.json', 'w', encoding='utf-8') as f:
                json.dump(json_entry, f, ensure_ascii=False, indent=4)

def compute_self_diffusion_coefficients_and_write_to_json(water_name,
                                                          molecule_names,
                                                          chi_water, 
                                                          molecule_abbreviations, 
                                                          water_abbreviation, 
                                                          mdanalysis_universe, 
                                                          time_between_frames,
                                                          res_rawdir,
                                                          name,
                                                          replica_number,
                                                          resdir,
                                                          temperature,
                                                          pressure,
                                                          total_number_molecules,
                                                          folder,
                                                          partsdir):
    """Computes self-diffusion coefficients and writes it to part JSON files.

    Args:
        water_name (str): The water name.
        molecule_names (list):  A list with all molecule names of the mixture.
        chi_water (float): Water mole fraction of the box.
        molecule_abbreviations (list): A list with all molecule abbreviations of the mixture.
        water_abbreviation (str): The water abbreviation that is used inside the PDB file.
        mdanalysis_universe (obj): MDAnalysis universe.
        time_between_frames (float): The time between frames in the trajectory.
        res_rawdir (str): Directory with raw, unprocessed results.
        name (str): System name.
        replica_number (int): Replica number.
        resdir (str): Path to the results directory.
        temperature (float): Temperature in kelvin at which the box was simulated.
        pressure (float): Pressure in bar at which the box was simulated.
        total_number_molecules (int): Total number of molecules inside the simulation box.
        folder (str): Path to replica folder.
        partsdir (str): Path to the parts directory.
    """

    if f'{water_name}' not in molecule_names:
            molecule_names.append(f'{water_name}')
    if chi_water != 0:
        molecule_abbreviations.append(f'{water_abbreviation}')
    if chi_water == 1:
        molecule_abbreviations = [water_abbreviation]
    for abbreviation in molecule_abbreviations:
        msd_time_axis, msd_timeseries = compute_msd_and_write_to_csv(mdanalysis_universe, abbreviation, time_between_frames, res_rawdir, name, replica_number)
        fit_slope = linear_fit_self_diffusion_coefficients_and_plot_msd(msd_time_axis, msd_timeseries, res_rawdir, name, abbreviation)
        write_self_diffusion_coefficients_to_json(resdir, 
                                                  temperature, 
                                                  pressure, 
                                                  chi_water, 
                                                  total_number_molecules, 
                                                  fit_slope, 
                                                  folder, 
                                                  replica_number, 
                                                  molecule_names,
                                                  molecule_abbreviations,
                                                  abbreviation,
                                                  partsdir,
                                                  name)
        del msd_time_axis
        del msd_timeseries

# This is the part that runs
mdanalysis_universe = load_mdanalysis_universe(f'{folder}/nvt_equilibration_end.pdb', f'{folder}/nve_production_unwrapped.dcd')
if mdanalysis_universe is not None:
    compute_self_diffusion_coefficients_and_write_to_json(water_name,
                                                          molecule_names,
                                                          chi_water, 
                                                          molecule_abbreviations, 
                                                          water_abbreviation, 
                                                          mdanalysis_universe, 
                                                          time_between_frames,
                                                          res_rawdir,
                                                          name,
                                                          replica_number,
                                                          resdir,
                                                          temperature,
                                                          pressure,
                                                          total_number_molecules,
                                                          folder,
                                                          partsdir)
else:
    with open(f'{resdir}/msd_fit_error_report.txt', 'a+') as f:
        f.write(f'{folder} could not be loaded as MDAnalysis universe.')