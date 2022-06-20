"""Module with useful functions for preparing a project.

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

def max_permitted_number_of_molecules(molecule_atom_numbers=[], molecule_molar_ratios=[]):
    """Calculates the maximum number of molecules that can fit inside a PDB file with the given parameters.

    Args:
        molecule_atom_numbers (list, optional): The number of atoms that the molecule consists of. Defaults to [].
        molecule_molar_ratios (list, optional): The proportion of the respective molecule type in the mixture. Defaults to [].

    Returns:
        int: The maximum number of molecules that can fit inside a PDB file with the given parameters.
    """

    PDB_MAX_ATOMS = 99999
    smallest_package_number_atoms = sum([n*r for n, r in zip(molecule_atom_numbers, molecule_molar_ratios)])
    max_smallest_package = PDB_MAX_ATOMS / smallest_package_number_atoms
    max_permitted_number = max_smallest_package * sum(molecule_molar_ratios)
    return int(max_permitted_number)

def estimate_init_box_side_length(total_number_molecules, molecule_molar_mass, density_kg_m3):
    """Estimates the initial box size for a molecule type.

    Args:
        total_number_molecules (int): Number of molecules in the planned box.
        molecule_molar_mass (float): Molecular mass of the molecule type in mol per gram.
        density_kg_m3 (int): Density of the molecule type in kilogram per cubic meter.

    Returns:
        float: Initial box size for a molecule type.
    """
    
    AVOGADRO_NUMBER = 6.022 * (10 ** 23)
    mass_kg = (molecule_molar_mass / AVOGADRO_NUMBER) * total_number_molecules * (10 ** (-3))
    volume_m3 = mass_kg / density_kg_m3
    side_lenght_m = volume_m3 ** (1 / 3)
    side_lenght_angstrom = side_lenght_m * (10 ** 10)
    return side_lenght_angstrom