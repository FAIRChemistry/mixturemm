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
import pandas as pd
import json
import statistics

# get arguments from the command-line this file was started with
workspace = sys.argv[1]
name = sys.argv[2]
half_npt_equilibration_csv_columns = int(sys.argv[3])
folder = f'{workspace}/{name}/npt_equilibration'
partsdir = f'{workspace}/results/parts'
total_number_molecules = int(sys.argv[4])
temperature = float(sys.argv[5])
pressure = float(sys.argv[6])
chi_water = float(sys.argv[7])
molecule_names_ = sys.argv[8]
molecule_names = molecule_names_.split(',')

def get_avg_density_and_volume(path_to_csv_file, half_npt_equilibration_csv_columns):
    """Extracts the average density and volume of a system from the csv file obtained through the StateDataReporter of OpenMM.

    Args:
        path_to_csv_file (str): Path to the csv file obtained by the StateDataReporter of OpenMM.
        half_npt_equilibration_csv_columns (int): Number of columns of half of the NpT equilibration csv file.

    Returns:
        float: Average density and volume.
    """

    npt_equilibration_df = pd.read_csv(f'{path_to_csv_file}')
    half_npt_equilibration_df = npt_equilibration_df.tail(n=half_npt_equilibration_csv_columns)
    density_npt_equilibration = half_npt_equilibration_df.loc[:, 'Density (g/mL)']
    volume_npt_equilibration = half_npt_equilibration_df.loc[:, 'Box Volume (nm^3)']
    average_density_npt_equilibration = statistics.mean(density_npt_equilibration)
    average_volume_npt_equilibration = statistics.mean(volume_npt_equilibration)
    stdev_density = statistics.stdev(density_npt_equilibration)
    stdev_volume = statistics.stdev(volume_npt_equilibration)
    return average_density_npt_equilibration, stdev_density, average_volume_npt_equilibration, stdev_volume

def write_density_json_part(average_density, stdev_density, average_volume, stdev_volume, temperature, pressure, chi_water, molecule_names, partsdir, name):
    """Writes density part JSON files containing the average density and volume of one system.

    Args:
        average_density (float): Average density of a system.
        average_volume (float): Average volume of a system.
        temperature (float): Temperature in kelvin at which the box was simulated.
        pressure (_float): Pressure in bar at which the box was simulated.
        chi_water (float): Water mole fraction of the box.
        molecule_names (list): List with the names of all mixture molecules.
        partsdir (str): Path to the parts directory.
        name (str): System name string.
    """

    json_entry = {
        'total_number_molecules' : total_number_molecules,
        'temperature' : temperature,
        'pressure' : pressure,
        'chi_water' : chi_water,
        'molecules' : molecule_names,
        'density' : average_density,
        'stdev_density' : stdev_density,
        'volume' : average_volume,
        'stdev_volume' : stdev_volume
    }
    with open(f'{partsdir}/{name}_density.json', 'w', encoding='utf-8') as f:
        json.dump(json_entry, f, ensure_ascii=False, indent=4)

# This is the part that runs
average_density, stdev_density, average_volume, stdev_volume = get_avg_density_and_volume(f'{folder}/npt_equilibration.csv', half_npt_equilibration_csv_columns)
write_density_json_part(average_density, stdev_density, average_volume, stdev_volume, temperature, pressure, chi_water, molecule_names, partsdir, name)