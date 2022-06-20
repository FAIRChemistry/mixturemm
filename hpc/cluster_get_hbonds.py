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

import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as hba
import json
import sys

# get arguments from the command-line this file was started with
workspace = sys.argv[1]
name = sys.argv[2]
replica_number = int(sys.argv[3])
folder = f'{workspace}/{name}/replica_{replica_number}'
resdir = f'{workspace}/results'
res_rawdir = f'{workspace}/results/raw_data'
partsdir = f'{workspace}/results/parts'
chi_water = float(sys.argv[4])
temperature = float(sys.argv[5])
pressure = float(sys.argv[6])
total_number_molecules = int(sys.argv[7])
molecule_names_ = sys.argv[8]
molecule_names = molecule_names_.split(',')
donors = sys.argv[9]
hydrogens = sys.argv[10]
acceptors = sys.argv[11]

# donors, hydrogens und acceptors as list ['resname HOH and name O', 'resname Me and name O']

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

def compute_hydrogen_bonds_and_write_to_json(mdanalysis_universe, 
                                             donors, 
                                             hydrogens, 
                                             acceptors,
                                             total_number_molecules,
                                             replica_number,
                                             temperature,
                                             pressure,
                                             chi_water,
                                             molecule_names,
                                             partsdir,
                                             name):
    """Computes the hydrogen bonds and writes them to part JSON files.

    Args:
        mdanalysis_universe (obj): MDAnalysis universe.
        donors (str): All hydrogen bond donors in MDAnalysis selection syntax.
        hydrogens (str): All hydrogen bond hydrogens in MDAnalysis selection syntax.
        acceptors (str): All hydrogen bond acceptors in MDAnalysis selection syntax.
        total_number_molecules (int): Total number of molecules inside the simulation box.
        replica_number (int): Replica number.
        temperature (float): Temperature in kelvin at which the box was simulated.
        pressure (float): Pressure in bar at which the box was simulated.
        chi_water (float): Water mole fraction of the box.
        molecule_names (list):  A list with all molecule names of the mixture.
        partsdir (str): Path to the parts directory.
        name (str): System name.
    """

    Hbonds = hba(universe=mdanalysis_universe,
                 donors_sel=donors,
                 hydrogens_sel=hydrogens,
                 acceptors_sel=acceptors,
                 d_a_cutoff=3.5,
                 d_h_a_angle_cutoff=150,
                 update_selections=False)
    Hbonds.run()
    try:
        hbonds_list = []
        for donor, acceptor, count in Hbonds.count_by_type():
            donor_resname, donor_type = donor.split(":")
            n_donors = mdanalysis_universe.select_atoms(f"resname {donor_resname} and type {donor_type}").n_atoms
            mean_count = 2 * int(count) / (Hbonds.n_frames * n_donors)
            hbonds_list.append(f"{donor} to {acceptor}: {mean_count:.2f}")
        json_entry = {
            'total_number_molecules' : total_number_molecules,
            'replica_number' : replica_number,
            'temperature' : temperature,
            'pressure' : pressure,
            'chi_water' : chi_water,
            'molecules' : molecule_names,
            'hydrogen_bonds' : hbonds_list
        }
        with open(f'{partsdir}/{name}_replica_{replica_number}_hbonds.json', 'w', encoding='utf-8') as f:
            json.dump(json_entry, f, ensure_ascii=False, indent=4)
    except:
        with open(f'{resdir}/hbonds_error_report.txt', 'a+') as f:
            f.write(f'{folder} could not be processed to hydrogen bonds per molecule.')
    del Hbonds

# This is the part that runs
mdanalysis_universe = load_mdanalysis_universe(f'{folder}/nvt_equilibration_end.pdb', f'{folder}/nve_production_wrapped.dcd')
if mdanalysis_universe is not None:
    compute_hydrogen_bonds_and_write_to_json(mdanalysis_universe, 
                                             donors, 
                                             hydrogens, 
                                             acceptors,
                                             total_number_molecules,
                                             replica_number,
                                             temperature,
                                             pressure,
                                             chi_water,
                                             molecule_names,
                                             partsdir,
                                             name)
else:
    with open(f'{resdir}/hbonds_error_report.txt', 'a+') as f:
        f.write(f'{folder} could not be loaded as MDAnalysis universe.')