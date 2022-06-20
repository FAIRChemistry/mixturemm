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

import json
import sys
import itertools
import matplotlib.pyplot as plt
import statistics
plt.switch_backend('agg')
plt.ioff()

# get arguments from the command-line this file was started with
chi_water_s__ = sys.argv[1]
chi_water_s_ = chi_water_s__.split(',')
chi_water_s = [float(f) for f in chi_water_s_]
temperature_s__ = sys.argv[2]
temperature_s_ = temperature_s__.split(',')
temperature_s = [float(f) for f in temperature_s_]
npt_equilibration_pressure_s__ = sys.argv[3]
npt_equilibration_pressure_s_ = npt_equilibration_pressure_s__.split(',')
npt_equilibration_pressure_s = [float(f) for f in npt_equilibration_pressure_s_]
workspace = sys.argv[4]
name_list_ = sys.argv[5]
molecule_names_ = sys.argv[6]
molecule_names = molecule_names_.split(',')
water_name = sys.argv[7]
name_list = name_list_.split(',')
description = sys.argv[8]
resdir = f'{workspace}/results'
partsdir = f'{workspace}/results/parts'

def write_density_json_header(description, resdir):
    """Writes the header for the density JSON file.

    Args:
        description (str): Path to the project description JSON file.
        resdir (str): Path to the results directory.
    """

    with open(f'{description}', 'r+', encoding='utf-8') as f:
        json_data = json.load(f)
        json_data['units'] = [
        {'temperature' : 'kelvin'},
        {'pressure' : 'bar'},
        {'density' : 'gram per cubic centimeter'},
        {'volume' : 'cubic nanometer'}
        ]
        json_data['densities'] = []
    with open(f'{resdir}/densities.json', 'w', encoding='utf-8') as f:
        json.dump(json_data, f, ensure_ascii=False, indent=4)

def scrape_and_add_density_parts(resdir, partsdir, name_list):
    """Reads and transfers all density part JSON files from the parts directory to the density JSON file.

    Args:
        resdir (str): Path to the results directory.
        partsdir (str): Path to the parts directory.
        name_list (list): List with all system names.
    """

    with open(f'{resdir}/densities.json', 'r+', encoding='utf-8') as f:
        json_data = json.load(f)
        for name in name_list:
            with open(f'{partsdir}/{name}_density.json', 'r', encoding='utf-8') as file:
                json_entry = json.load(file)
            json_data['densities'].append(json_entry)
        f.seek(0)
        json.dump(json_data, f, ensure_ascii=False, indent=4)

def write_averaged_density_over_boxes_to_json(resdir, temperature_s, npt_equilibration_pressure_s, chi_water_s):
    """Computes the average density over all box sizes and writes it to the density JSON file.

    Args:
        resdir (str): Path to the results directory.
        temperature_s (list): A list with all different temperatures in kelvin at which the boxes were simulated.
        npt_equilibration_pressure_s (list): A list with all pressures in bar at which the boxes were simulated.
        chi_water_s (list): A list with all different molar fractions of water that were simulated.
    """

    with open(f'{resdir}/densities.json', 'r+', encoding='utf-8') as f:
        json_data = json.load(f)
        json_data['averaged densities'] = []
        parameter_list = [temperature_s, npt_equilibration_pressure_s, chi_water_s]
        combination_list = list(itertools.product(*parameter_list))
        for combination in combination_list:
            density_list = []
            density_stdev_list = []
            temperature, pressure, chi_water = combination
            for entry in json_data['densities']:
                if combination == (entry['temperature'], entry['pressure'], entry['chi_water']):
                    density_list.append(entry['density'])
                    density_stdev_list.append(entry['stdev_density'])
            density_mean = statistics.mean(density_list)
            density_stdev = statistics.mean(density_stdev_list)
            json_entry = {
                'temperature' : temperature,
                'pressure' : pressure,
                'chi_water' : chi_water,
                'molecules' : molecule_names,
                'average_density' : density_mean,
                'stdev_density' : density_stdev
            }
            json_data['averaged densities'].append(json_entry)
        f.seek(0)
        json.dump(json_data, f, ensure_ascii=False, indent=4)

def plot_average_densities_over_chi_water(resdir, npt_equilibration_pressure_s, temperature_s, molecule_names):
    """Plots the average densities over the water mole fractions.

    Args:
        resdir (str): Path to the results directory.
        npt_equilibration_pressure_s (list): A list with all pressures in bar at which the boxes were simulated.
        temperature_s (list): A list with all different temperatures in kelvin at which the boxes were simulated.
        molecule_names (list): A list with all molecule names of the mixture.
    """

    with open(f'{resdir}/densities.json', 'r', encoding='utf-8') as f:
        json_data = json.load(f)
    for pressure in npt_equilibration_pressure_s:
        for temperature in temperature_s:
            x_values = []
            y_values = []
            y_error = []
            for entry in json_data['averaged densities']:
                if temperature == entry['temperature'] and pressure == entry['pressure']:
                    x_values.append(entry['chi_water'])
                    y_values.append(entry['average_density'])
                    y_error.append(entry['stdev_density'])
            plt.errorbar(x_values, y_values, yerr=y_error, label=f'{temperature} K')
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
        title = molecule_names
        if f'{water_name}' not in molecule_names:
            molecule_names.append(f'{water_name}')
            title_molecule_list = molecule_names
            title = '_'.join(title_molecule_list)
            title_molecule_list.pop()
        plt.title(f'Density {title} at {pressure} bar')
        plt.xlabel(r'$\chi_{Water}$')
        plt.ylabel(r'density in $\frac{g}{cm^3}$')
        plt.savefig(f'{resdir}/Density_{title}_at_{pressure}_bar.png', dpi=300, bbox_inches='tight')
        plt.clf()

# This is the part that runs
write_density_json_header(description, resdir)
scrape_and_add_density_parts(resdir, partsdir, name_list)
write_averaged_density_over_boxes_to_json(resdir, temperature_s, npt_equilibration_pressure_s, chi_water_s)
if len(chi_water_s) > 2:
    plot_average_densities_over_chi_water(resdir, npt_equilibration_pressure_s, temperature_s, molecule_names)