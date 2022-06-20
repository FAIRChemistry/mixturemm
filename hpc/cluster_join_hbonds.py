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
import itertools
import matplotlib.pyplot as plt
import sys
import re
import statistics
plt.switch_backend('agg')
plt.ioff()

workspace = sys.argv[1]
name_list_ = sys.argv[2]
number_dict__ = sys.argv[3]
number_dict_ = dict(x.split('=') for x in number_dict__.split(','))
number_dict = dict([tot, rep] for tot, rep in number_dict_.items())
name_list = name_list_.split(',')
resdir = f'{workspace}/results'
res_rawdir = f'{workspace}/results/raw_data'
partsdir = f'{workspace}/results/parts'
description = sys.argv[4]
chi_water_s__ = sys.argv[5]
chi_water_s_ = chi_water_s__.split(',')
chi_water_s = [float(f) for f in chi_water_s_]
temperature_s__ = sys.argv[6]
temperature_s_ = temperature_s__.split(',')
temperature_s = [float(f) for f in temperature_s_]
npt_equilibration_pressure_s__ = sys.argv[7]
npt_equilibration_pressure_s_ = npt_equilibration_pressure_s__.split(',')
npt_equilibration_pressure_s = [float(f) for f in npt_equilibration_pressure_s_]
molecule_names_ = sys.argv[8]
molecule_names = molecule_names_.split(',')
water_name = sys.argv[9]
molecule_abbreviations_ = sys.argv[10]
molecule_abbreviations = molecule_abbreviations_.split(',')
water_abbreviation = sys.argv[11]
replica_start = int(sys.argv[12])

def flatten_list(list_to_be_flattend):
    """Flattens a list of lists.

    Args:
        list_to_be_flattend (list): List that should be flattend.

    Returns:
        list: Flattend list.
    """

    return list(itertools.chain.from_iterable(list_to_be_flattend))

with open(f'{description}', 'r+', encoding='utf-8') as f:
    json_data = json.load(f)
    json_data['units'] = [
    {'temperature' : 'kelvin'},
    {'pressure' : 'bar'},
    {'hydrogen_bonds' : 'count per molecule'}
    ]
    json_data['hydrogen bonds'] = []
with open(f'{resdir}/hydrogen_bonds.json', 'w', encoding='utf-8') as f:
    json.dump(json_data, f, ensure_ascii=False, indent=4)
with open(f'{resdir}/hydrogen_bonds.json', 'r+', encoding='utf-8') as f:
    json_data = json.load(f)
    combination_list_ = []
    for total_num, rep_count in number_dict.items():
        new_rep_numbers = list(range(replica_start, (replica_start + int(rep_count))))
        filtered_names = [entry for entry in name_list if entry.split('-')[0] == str(total_num)]
        parameter_list = [filtered_names, new_rep_numbers]
        combination_list_part = list(itertools.product(*parameter_list))
        combination_list_.append(combination_list_part)
    combination_list = flatten_list(combination_list_)
    for combination in combination_list:
        name, number = combination
        try:
            with open(f'{partsdir}/{name}_replica_{number}_hbonds.json', 'r', encoding='utf-8') as file:
                json_entry = json.load(file)
            json_data['hydrogen bonds'].append(json_entry)
        except:
            pass
    f.seek(0)
    json.dump(json_data, f, ensure_ascii=False, indent=4)
with open(f'{resdir}/hydrogen_bonds.json', 'r+', encoding='utf-8') as f:
    json_data = json.load(f)
    for pressure in npt_equilibration_pressure_s:
        for temperature in temperature_s:
            x_values = []
            y_values = []
            zero_list = []
            one_list = []
            for entry in json_data['hydrogen bonds']:
                if temperature == entry['temperature'] and entry['chi_water'] == 0 and pressure == entry['pressure']:
                    zero_list.append(entry['hydrogen_bonds'])
                if temperature == entry['temperature'] and entry['chi_water'] == 1 and pressure == entry['pressure']:
                    one_list.append(entry['hydrogen_bonds'])
            zero_sum_list = []
            one_sum_list = []
            for element in zero_list:
                joined_zero_string = '-'.join(element)
                zero_sum_list.append(re.findall("\d+\.\d+", joined_zero_string))
            for element in one_list:
                joined_one_string = '-'.join(element)
                one_sum_list.append(re.findall("\d+\.\d+", joined_one_string))
            zero_sum_list_flat = [item for sublist in zero_sum_list for item in sublist]
            one_sum_list_flat = [item for sublist in one_sum_list for item in sublist]
            zero_sum_floats = [float(item) for item in zero_sum_list_flat]
            one_sum_floats = [float(item) for item in one_sum_list_flat]
            zero_value = statistics.mean(zero_sum_floats)
            one_value = statistics.mean(one_sum_floats)
            for chi_value in chi_water_s:
                chi_list = []
                for entry in json_data['hydrogen bonds']:
                    if temperature == entry['temperature'] and pressure == entry['pressure'] and chi_value == entry["chi_water"]:
                        joined_hbonds_string = '-'.join(entry['hydrogen_bonds'])
                        chi_list.append(re.findall("\d+\.\d+", joined_hbonds_string))
                        chi_floats = [list(map(float, item)) for item in chi_list]
                        chi_list_flat = [sum(sublist) for sublist in chi_floats]
                        x_values.append(chi_value)
                        y_values.append(statistics.mean(chi_list_flat) - (chi_value * (one_value - zero_value) + zero_value))
            plt.plot(x_values, y_values, linestyle='-', marker='o', label=f'{temperature} K')
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
        title = molecule_names
        if f'{water_name}' not in molecule_names:
            molecule_names.append(f'{water_name}')
            title_molecule_list = molecule_names
            title = '_'.join(title_molecule_list)
            title_molecule_list.pop()
        plt.title(f'Excess hydrogen bonds {title} at {pressure} bar')
        plt.xlabel(r'$\chi_{Water}$')
        plt.ylabel(r'excess density')
        plt.savefig(f'{resdir}/Excess_hydrogen_bonds_{title}_at_{pressure}_bar.png', dpi=300, bbox_inches='tight')
        plt.clf()
with open(f'{resdir}/hydrogen_bonds.json', 'r+', encoding='utf-8') as f:
    json_data = json.load(f)
    for pressure in npt_equilibration_pressure_s:
        for temperature in temperature_s:
            for ab in molecule_abbreviations:
                x_values = []
                y_values = []
                filtered_chi_values = [chi for chi in chi_water_s if chi != 1]
                for chi_value in filtered_chi_values:
                    chi_list = []
                    for entry in json_data['hydrogen bonds']:
                        if temperature == entry['temperature'] and pressure == entry['pressure'] and chi_value == entry["chi_water"]:
                            filtered_hbonds_list = [item for item in entry['hydrogen_bonds'] if item.startswith(f'{ab}:')]
                            joined_hbonds_string = '-'.join(filtered_hbonds_list)
                            chi_list.append(re.findall("\d+\.\d+", joined_hbonds_string))
                            chi_floats = [list(map(float, item)) for item in chi_list]
                            chi_list_flat = [sum(sublist) for sublist in chi_floats]
                            x_values.append(chi_value)
                            y_values.append(statistics.mean(chi_list_flat))
                plt.plot(x_values, y_values, linestyle='-', marker='o', label=f'{temperature} K {ab}')
            x_values = []
            y_values = []
            filtered_chi_values = [chi for chi in chi_water_s if chi != 0]
            for chi_value in filtered_chi_values:
                chi_list = []
                for entry in json_data['hydrogen bonds']:
                    if temperature == entry['temperature'] and pressure == entry['pressure'] and chi_value == entry["chi_water"]:
                        filtered_hbonds_list = [item for item in entry['hydrogen_bonds'] if item.startswith(f'{water_abbreviation}:')]
                        joined_hbonds_string = '-'.join(filtered_hbonds_list)
                        chi_list.append(re.findall("\d+\.\d+", joined_hbonds_string))
                        chi_floats = [list(map(float, item)) for item in chi_list]
                        chi_list_flat = [sum(sublist) for sublist in chi_floats]
                        x_values.append(chi_value)
                        y_values.append(statistics.mean(chi_list_flat))
            plt.plot(x_values, y_values, linestyle='-', marker='o', label=f'{temperature} K {water_abbreviation}')
        plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
        title_molecule_list = molecule_names.append(water_name)
        title = '_'.join(title_molecule_list)
        title_molecule_list.pop()
        plt.title(f'Hydrogen bonds per molecule {title} at {pressure} bar')
        plt.xlabel(r'$\chi_{Water}$')
        plt.ylabel(r'hydrogen bonds per molecule')
        plt.savefig(f'{resdir}/Hydrogen_bonds_per_molecule_{title}_at_{pressure}_bar.png', dpi=300, bbox_inches='tight')
        plt.clf()