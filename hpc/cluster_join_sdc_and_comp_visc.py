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

import numpy as np
from scipy.stats import linregress
import json
import itertools
import matplotlib.pyplot as plt
import math
import statistics
import sys
plt.switch_backend('agg')
plt.ioff()

# get arguments from the command-line this file was started with
K_B = 1.380649*(10**(-23))
XI = 2.837297
PI = math.pi
R = 8.31446261815324
COVERSION_FACTOR_PAS = 10 ** 9
total_number_molecules__ = sys.argv[1]
total_number_molecules_ = total_number_molecules__.split(',')
total_number_molecules = [int(f) for f in total_number_molecules_]
chi_water_s__ = sys.argv[2]
chi_water_s_ = chi_water_s__.split(',')
chi_water_s = [float(f) for f in chi_water_s_]
temperature_s__ = sys.argv[3]
temperature_s_ = temperature_s__.split(',')
temperature_s = [float(f) for f in temperature_s_]
npt_equilibration_pressure_s__ = sys.argv[4]
npt_equilibration_pressure_s_ = npt_equilibration_pressure_s__.split(',')
npt_equilibration_pressure_s = [float(f) for f in npt_equilibration_pressure_s_]
workspace = sys.argv[5]
name_list_ = sys.argv[6]
molecule_names_ = sys.argv[7]
molecule_names = molecule_names_.split(',')
water_name = sys.argv[8]
number_dict__ = sys.argv[9]
number_dict_ = dict(x.split('=') for x in number_dict__.split(','))
number_dict = dict([tot, rep] for tot, rep in number_dict_.items())
name_list = name_list_.split(',')
resdir = f'{workspace}/results'
res_rawdir = f'{workspace}/results/raw_data'
partsdir = f'{workspace}/results/parts'
molecule_abbreviations_ = sys.argv[10]
molecule_abbreviations = molecule_abbreviations_.split(',')
water_abbreviation = sys.argv[11]
description = sys.argv[12]
replica_start = int(sys.argv[13])

def flatten_list(list_to_be_flattend):
    """Flattens a list of lists.

    Args:
        list_to_be_flattend (list): List that should be flattend.

    Returns:
        list: Flattend list.
    """

    return list(itertools.chain.from_iterable(list_to_be_flattend))

def box_size_sort_key(list_element):
    return int(list_element[0].split('-')[0])

def write_sdc_json_header(description, resdir):
    """Writes the header for the self-diffusion coefficients JSON file.

    Args:
        description (str): Path to the project description JSON file.
        resdir (str): Path to the results directory.
    """

    with open(f'{description}', 'r+', encoding='utf-8') as f:
        json_data = json.load(f)
        json_data['units'] = [
        {'temperature' : 'kelvin'},
        {'pressure' : 'bar'},
        {'self_diffusion_coefficient' : 'square meter per second'}
        ]
        json_data['self diffusion coefficients'] = []
    with open(f'{resdir}/self_diffusion_coefficients.json', 'w', encoding='utf-8') as f:
        json.dump(json_data, f, ensure_ascii=False, indent=4)

def scrape_and_add_sdc_parts(resdir, 
                             water_abbreviation, 
                             molecule_abbreviations, 
                             number_dict, 
                             replica_start, 
                             partsdir):
    """Reads and transfers all self-diffusion coefficient part JSON files from the parts directory to the self-diffusion coefficients JSON file.

    Args:
        resdir (str): Path to the results directory.
        water_abbreviation (str): The water abbreviation that is used inside the PDB file.
        molecule_abbreviations (list): A list with all molecule abbreviations of the mixture.
        number_dict (dict): Contains all replica amounts of the systems specified by total number of molecules.
        replica_start (int): The number from which on the replica numbering starts.
        partsdir (str): Path to the parts directory.
    """

    with open(f'{resdir}/self_diffusion_coefficients.json', 'r+', encoding='utf-8') as f:
        json_data = json.load(f)
        if f'{water_abbreviation}' not in molecule_abbreviations:
            molecule_abbreviations.append(f'{water_abbreviation}')
        combination_list___ = []
        for total_num, rep_count in number_dict.items():
            new_rep_numbers = list(range(replica_start, (replica_start + int(rep_count))))
            filtered_names = [entry for entry in name_list if entry.split('-')[0] == str(total_num)]
            parameter_list = [filtered_names, new_rep_numbers, molecule_abbreviations]
            combination_list_part = list(itertools.product(*parameter_list))
            combination_list___.append(combination_list_part)
        combination_list__ = flatten_list(combination_list___)
        combination_list_ = [t for t in combination_list__ if not (t[0].split('-')[2] == '0' and t[2] == 'HOH')]
        combination_list = [t for t in combination_list_ if not (t[0].split('-')[2] == '1' and t[2] != 'HOH')]
        combination_list.sort(key=box_size_sort_key)
        for combination in combination_list:
            name, number, abb = combination
            try:
                with open(f'{partsdir}/{name}_{abb}_replica_{number}_sdc.json', 'r', encoding='utf-8') as file:
                    json_entry = json.load(file)
                json_data['self diffusion coefficients'].append(json_entry)
            except:
                pass
        f.seek(0)
        json.dump(json_data, f, ensure_ascii=False, indent=4)

def write_averaged_sdc_over_boxes_to_json(resdir, 
                                          water_abbreviation, 
                                          molecule_abbreviations, 
                                          temperature_s, 
                                          npt_equilibration_pressure_s, 
                                          chi_water_s, 
                                          total_number_molecules,
                                          molecule_names):
    """Computes the average self-diffusion coefficient over all replicas of a system and writes it to the self-diffusion coefficients JSON file.

    Args:
        resdir (str): Path to the results directory.
        water_abbreviation (str): The water abbreviation that is used inside the PDB file.
        molecule_abbreviations (list): A list with all molecule abbreviations of the mixture.
        temperature_s (list): A list with all different temperatures in kelvin at which the boxes were simulated.
        npt_equilibration_pressure_s (list): A list with all pressures in bar at which the boxes were simulated.
        chi_water_s (list): A list with all different molar fractions of water that were simulated.
        total_number_molecules (list): A list with the number of molecules that were placed in the simulation box.
        molecule_names (list): A list with all molecule names of the mixture.
    """

    with open(f'{resdir}/self_diffusion_coefficients.json', 'r+', encoding='utf-8') as f:
        json_data = json.load(f)
        json_data['averaged self diffusion coefficients'] = []
        if f'{water_abbreviation}' not in molecule_abbreviations:
            molecule_abbreviations.append(f'{water_abbreviation}')
        parameter_list = [temperature_s, npt_equilibration_pressure_s, chi_water_s, molecule_abbreviations, total_number_molecules]
        combination_list__ = list(itertools.product(*parameter_list))
        combination_list_ = [t for t in combination_list__ if not (t[2] == 0 and t[3] == 'HOH')]
        combination_list = [t for t in combination_list_ if not (t[2] == 1 and t[3] != 'HOH')]
        for combination in combination_list:
            self_diffusion_coefficient_list = []
            temperature, pressure, chi_water , abbreviation, number = combination
            for entry in json_data['self diffusion coefficients']:
                if combination == (entry['temperature'], entry['pressure'], entry['chi_water'], entry['abbreviation'], entry['total_number_molecules']):
                    self_diffusion_coefficient_list.append(entry['self_diffusion_coefficient'])
                    volume = entry['volume']
            self_diffusion_coefficient_mean = statistics.mean(self_diffusion_coefficient_list)
            if len(self_diffusion_coefficient_list) == 1:
                self_diffusion_coefficient_stdev = 0
            else:
                self_diffusion_coefficient_stdev = statistics.stdev(self_diffusion_coefficient_list)
            json_entry = {
                'total_number_molecules' : number,
                'temperature' : temperature,
                'pressure' : pressure,
                'chi_water' : chi_water,
                'volume' : volume,
                'molecules' : molecule_names,
                'abbreviations' : molecule_abbreviations,
                'abbreviation' : abbreviation,
                'average_self_diffusion_coefficient' : self_diffusion_coefficient_mean,
                'stdev_self_diffusion_coefficient' : self_diffusion_coefficient_stdev
            }
            json_data['averaged self diffusion coefficients'].append(json_entry)
        f.seek(0)
        json.dump(json_data, f, ensure_ascii=False, indent=4)

def write_visc_json_header(description, resdir):
    """Writes the header for the viscosities JSON file.

    Args:
        description (str): Path to the project description JSON file.
        resdir (str): Path to the results directory.
    """

    with open(f'{description}', 'r+', encoding='utf-8') as f:
        json_data = json.load(f)
        json_data['units'] = [
        {'temperature' : 'kelvin'},
        {'pressure' : 'bar'},
        {'self_diffusion_coefficient' : 'meter per square seconds'},
        {'viscosity' : 'pascal seconds'},
        {'e eta' : 'kilojoule per mole'}
        ]
        json_data['viscosities'] = []
    with open(f'{resdir}/viscosities.json', 'w', encoding='utf-8') as f:
        json.dump(json_data, f, ensure_ascii=False, indent=4)

def get_values_for_comp_and_plot_visc(data, combination):
    """Gets the self-diffusion coefficients needed for viscosity computation from the self-diffusion coefficients JSON.

    Args:
        data (dict): The data of the self-diffusion coefficients JSON.
        combination (tuple): Parameters of self-diffusion coefficient measurements for a system.

    Returns:
        list: self-diffusion coefficients over inverse box size for viscosity fitting.
    """

    x_values = []
    y_values = []
    for entry in data['self diffusion coefficients']:
        if combination == (entry['temperature'], entry['pressure'], entry['chi_water'], entry['abbreviation']):
            x_values.append(1 / ((entry['volume']) ** (1 / 3)))
            y_values.append(entry['self_diffusion_coefficient'])
    x_plot = []
    y_plot = []
    y_error_plot = []
    for entry in data['averaged self diffusion coefficients']:
        if combination == (entry['temperature'], entry['pressure'], entry['chi_water'], entry['abbreviation']):
            x_plot.append(1 / ((entry['volume']) ** (1 / 3)))
            y_plot.append(entry['average_self_diffusion_coefficient'])
            y_error_plot.append(entry['stdev_self_diffusion_coefficient'])
    return x_values, y_values, x_plot, y_plot, y_error_plot

def comp_visc_and_write_to_json(combination, 
                                x_values, 
                                y_values, 
                                K_B, 
                                XI, 
                                PI, 
                                COVERSION_FACTOR_PAS, 
                                resdir, 
                                molecule_names,
                                molecule_abbreviations):
    """Computes the viscosity and writes it to the viscosities JSON file.

    Args:
        combination (tuple): Parameters of self-diffusion coefficient measurements for a system.
        x_values (list): Inverse box sizes.
        y_values (list): Self-diffusion coefficients.
        K_B (float): Boltzmann constant.
        XI (float): Dimensionless constant for cubic boxes.
        PI (float): Pi.
        COVERSION_FACTOR_PAS (float): Conversion factor for viscosity in pascal seconds.
        resdir (str): Path to the results directory.
        molecule_names (list): A list with all molecule names of the mixture.
        molecule_abbreviations (list): A list with all molecule abbreviations of the mixture.

    Returns:
        float: Y-axis intercept and slope of the viscosity fit.
    """

    temperature, pressure, chi_water, abbreviation = combination
    linear_model = linregress(x_values, y_values)
    viscosity = ((- K_B * temperature * XI) / (linear_model.slope * 6 * PI)) * COVERSION_FACTOR_PAS
    viscosity_stderr = ((- K_B * temperature * XI) / (linear_model.stderr * 6 * PI)) * COVERSION_FACTOR_PAS
    self_diffusion_coefficient_infinite = linear_model.intercept
    self_diffusion_coefficient_infinite_stderr = linear_model.intercept_stderr
    with open(f'{resdir}/viscosities.json', 'r+', encoding='utf-8') as f:
        json_data = json.load(f)
        json_entry = {
            'temperature' : temperature,
            'pressure' : pressure,
            'chi_water' : chi_water,
            'molecules' : molecule_names,
            'abbreviations' : molecule_abbreviations,
            'abbreviation' : abbreviation,
            'viscosity' : viscosity,
            'viscosity_stderr' : viscosity_stderr,
            'self_diffusion_coefficient_infinite' : self_diffusion_coefficient_infinite,
            'self_diffusion_coefficient_infinite_stderr' : self_diffusion_coefficient_infinite_stderr
        }
        json_data['viscosities'].append(json_entry)
        f.seek(0)
        json.dump(json_data, f, ensure_ascii=False, indent=4)
    return linear_model.slope, linear_model.intercept

def plot_visc_comp(combination, molecule_names, water_name, x_plot, y_plot, y_error_plot, intercept, slope):
    """Plots the viscosity fit.

    Args:
        combination (tuple): Parameters of self-diffusion coefficient measurements for a system.
        molecule_names (list): A list with all molecule names of the mixture.
        water_name (str): The water name.
        x_plot (list): Inverse box sizes.
        y_plot (list): Average self-diffusion coefficients.
        y_error_plot (list): Standard deviation of average self-diffusion coefficients.
        intercept (float): Y-axis intercept of the viscosity fit.
        slope (float): Slope of the viscosity fit.
    """

    temperature, pressure, chi_water, abbreviation = combination
    title = molecule_names
    if f'{water_name}' not in molecule_names:
        molecule_names.append(f'{water_name}')
        title_molecule_list = molecule_names
        title = '_'.join(title_molecule_list)
        title_molecule_list.pop()
    plt.errorbar(x_plot, y_plot, yerr=y_error_plot, linestyle='None')
    x_line = np.arange(0, max(x_plot)+0.1, 0.001)
    plt.plot(x_line, intercept + slope * x_line, color='red')
    plt.title(r'$D_{\infty}$ fit')
    plt.xlabel(r'$\frac{1}{L}$')
    plt.ylabel(r'D in $\frac{m^2}{s}$')
    plt.savefig(f'{res_rawdir}/Visc_comp_{title}_{abbreviation}_at_{pressure}_bar.png', dpi=300)
    plt.clf()

def comp_all_visc(resdir, 
                  water_abbreviation, 
                  molecule_abbreviations,
                  temperature_s,
                  npt_equilibration_pressure_s,
                  chi_water_s,
                  K_B, 
                  XI, 
                  PI, 
                  COVERSION_FACTOR_PAS,
                  molecule_names,
                  water_name):
    """Computes the viscosity by linearly fitting the self-diffusion coefficients over their inverse box size.

    Args:
        resdir (str): Path to the results directory.
        water_abbreviation (str): The water abbreviation that is used inside the PDB file.
        molecule_abbreviations (list): A list with all molecule abbreviations of the mixture.
        temperature_s (list): A list with all different temperatures in kelvin at which the boxes were simulated.
        npt_equilibration_pressure_s (list): A list with all pressures in bar at which the boxes were simulated.
        chi_water_s (list): A list with all different molar fractions of water that were simulated.
        K_B (float): Boltzmann constant.
        XI (float): Dimensionless constant for cubic boxes.
        PI (float): Pi.
        COVERSION_FACTOR_PAS (float): Conversion factor for viscosity in pascal seconds.
        molecule_names (list): A list with all molecule names of the mixture.
        water_name (str): The water name.
    """

    with open(f'{resdir}/self_diffusion_coefficients.json', 'r', encoding='utf-8') as f:
        data = json.load(f)
    if f'{water_abbreviation}' not in molecule_abbreviations:
        molecule_abbreviations.append(f'{water_abbreviation}')
    parameter_list = [temperature_s, npt_equilibration_pressure_s, chi_water_s, molecule_abbreviations]
    combination_list__ = list(itertools.product(*parameter_list))
    combination_list_ = [t for t in combination_list__ if not (t[2] == 0 and t[3] == 'HOH')]
    combination_list = [t for t in combination_list_ if not (t[2] == 1 and t[3] != 'HOH')]
    for combination in combination_list:
        x_values, y_values, x_plot, y_plot, y_error_plot = get_values_for_comp_and_plot_visc(data, combination)
        slope, intercept = comp_visc_and_write_to_json(combination, 
                                                       x_values, 
                                                       y_values, 
                                                       K_B, 
                                                       XI, 
                                                       PI, 
                                                       COVERSION_FACTOR_PAS, 
                                                       resdir, 
                                                       molecule_names,
                                                       molecule_abbreviations)
        plot_visc_comp(combination, molecule_names, water_name, x_plot, y_plot, y_error_plot, intercept, slope)
        
def plot_viscosities_over_chi_water(resdir, 
                                    water_abbreviation, 
                                    molecule_abbreviations, 
                                    npt_equilibration_pressure_s,
                                    temperature_s,
                                    molecule_names,
                                    water_name):
    """Plots viscosities over their respective water mole fractions.

    Args:
        resdir (str): Path to the results directory.
        water_abbreviation (str): The water abbreviation that is used inside the PDB file.
        molecule_abbreviations (list): A list with all molecule abbreviations of the mixture.
        npt_equilibration_pressure_s (list): A list with all pressures in bar at which the boxes were simulated.
        temperature_s (list): A list with all different temperatures in kelvin at which the boxes were simulated.
        molecule_names (list): A list with all molecule names of the mixture.
        water_name (str): The water name.
    """

    with open(f'{resdir}/viscosities.json', 'r', encoding='utf-8') as f:
        json_data = json.load(f)
    if f'{water_abbreviation}' not in molecule_abbreviations:
        molecule_abbreviations.append(f'{water_abbreviation}')
    for pressure in npt_equilibration_pressure_s:
        for abbreviation in molecule_abbreviations:
            for temperature in temperature_s:
                x_values = []
                y_values = []
                y_error = []
                for entry in json_data['viscosities']:
                    if temperature == entry['temperature'] and abbreviation == entry['abbreviation'] and pressure == entry['pressure']:
                        x_values.append(entry['chi_water'])
                        y_values.append(entry['viscosity'])
                        y_error.append(entry['viscosity_stderr'])
                plt.errorbar(x_values, y_values, yerr=y_error, label=f'{temperature} K')
            plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
            title = molecule_names
            if f'{water_name}' not in molecule_names:
                molecule_names.append(f'{water_name}')
                title_molecule_list = molecule_names
                title = '_'.join(title_molecule_list)
                title_molecule_list.pop()
            plt.title(f'Viscosity {title} {abbreviation} at {pressure} bar')
            plt.xlabel(r'$\chi_{Water}$')
            plt.ylabel(r'viscosity in Pas')
            plt.savefig(f'{resdir}/Viscosity_{title}_{abbreviation}_at_{pressure}_bar.png', dpi=300, bbox_inches='tight')
            plt.clf()

# This is the part that runs
write_sdc_json_header(description, resdir)
scrape_and_add_sdc_parts(resdir, 
                         water_abbreviation, 
                         molecule_abbreviations, 
                         number_dict, 
                         replica_start, 
                         partsdir)
write_averaged_sdc_over_boxes_to_json(resdir, 
                                      water_abbreviation, 
                                      molecule_abbreviations, 
                                      temperature_s, 
                                      npt_equilibration_pressure_s, 
                                      chi_water_s, 
                                      total_number_molecules,
                                      molecule_names)
write_visc_json_header(description, resdir)
comp_all_visc(resdir, 
              water_abbreviation, 
              molecule_abbreviations,
              temperature_s,
              npt_equilibration_pressure_s,
              chi_water_s,
              K_B, 
              XI, 
              PI, 
              COVERSION_FACTOR_PAS,
              molecule_names,
              water_name)
if len(chi_water_s) > 2:
    plot_viscosities_over_chi_water(resdir, 
                                    water_abbreviation, 
                                    molecule_abbreviations, 
                                    npt_equilibration_pressure_s,
                                    temperature_s,
                                    molecule_names,
                                    water_name)