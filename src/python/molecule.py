import json
import numpy as np
import re

from basisfunction import BasisFunction
from atom import Atom
from orbital import Orbital

class Molecule:
    def __init__(self):
        self.atoms = list()
        self.basis_functions = list()
        self.orbitals = list()
        self.density_matrix = np.array([])

        self.json_data = None
        self.molden_data = None

        self.x_min = None
        self.x_max = None
        self.y_min = None
        self.y_max = None
        self.z_min = None
        self.z_max = None

    @classmethod
    def from_json(cls, json_file_name):
        instance = cls()

        json_file = open(json_file_name)
        instance.json_data = json.load(json_file)

        instance.read_basis_functions_from_json()
        instance.read_atoms_from_json()
        instance.read_densities_from_json()

        return instance

    @classmethod
    def from_molden(cls, molden_file_name):
        instance = cls()
        print("Loading a molden file ...")

        molden_file = open(molden_file_name)
        instance.molden_data = molden_file.read()
        print("Reading atom data ...")
        instance.read_atoms_from_molden()
        print("Reading GTO data ...")
        instance.read_basis_functions_from_molden()
        print("Reading MO data ...")
        instance.read_orbitals_from_molden()

        print("A molden file has been read!")

        # print(instance.molden_data)

        return instance
        
    def read_basis_functions_from_json(self):
        orbital_labels_list = self.json_data["Molecule"]["MolecularOrbitals"]["OrbitalLabels"]

        coefficients_list = list()
        exponents_list = list()
        shell_list = list()
        positions_list = list()

        for atom in self.json_data["Molecule"]["Atoms"]:
            # print("\tAtom number ", atom["Idx"])
            # print("Coefficients")

            for func in atom["Basis"]:
                # print(func["Coefficients"])
                coefficients_list.append(func["Coefficients"])

                # print("Exponents")
                # print(func["Exponents"])
                exponents_list.append(func["Exponents"])

                # print("Shell")
                # print(func["Shell"])
                shell_list.append(func["Shell"])
                
                positions_list.append(atom["Coords"])

        if not len(coefficients_list) == len(exponents_list) == len(shell_list) == len(positions_list):
            raise Exception("Class Molecule. The coefficients_list, exponents_list, shell_list and positions_list have different shapes")

        for coefficients, exponents, shell, position in zip(coefficients_list, exponents_list, shell_list, positions_list):
            if shell == "s":
                self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position))
            elif shell == "p":
                self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="z"))
                self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="x"))
                self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="y"))
            elif shell == "d":
                self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="z2"))
                self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="xz"))
                self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="yz"))
                self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="x2y2"))
                self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="xy"))

        if not len(orbital_labels_list) == len(self.basis_functions):
            raise Exception("Class Molecule. len(orbital_labels_list) is not equal len(self.basis_functions)")

        for orbital_label, basis_function in zip(orbital_labels_list, self.basis_functions):
            basis_function.label=orbital_label
            # print(basis_function)

    def read_atoms_from_json(self):
        for atom in self.json_data["Molecule"]["Atoms"]:
            self.atoms.append(Atom(atom["Coords"], atom["ElementLabel"], atom["Idx"], atom["NuclearCharge"]))
        self.update_molecule_limits()

    def read_densities_from_json(self):
        self.density_matrix = np.array(self.json_data["Molecule"]["Densities"]["scfp"])

        if not self.density_matrix.shape[0] == self.density_matrix.shape[1]:
            raise Exception("Class Molecule. Something is wrong with the scfp density")

        if not self.density_matrix.shape[0] == len(self.basis_functions):
            raise Exception("Class Molecule. The number of basis functions is inconsistent with the size of the density matrix")

        # print(np.array2string(self.density_matrix, separator=',', max_line_width=np.inf))
        # print(np.array_equal(self.density_matrix, self.density_matrix.T))

    def read_atoms_from_molden(self):
        # finding the [Atoms] block in the file
        start_block = "[Atoms] AU"
        end_block = "[GTO]"
        start_index = self.molden_data.find(start_block) + len(start_block)
        if start_index == -1:
            raise Exception("Class Molecule. There is no [Atoms] block in the file")

        end_index = self.molden_data.find(end_block, start_index)
        if end_index == -1:
            fragment = self.molden_data[start_index:]
        else:
            fragment = self.molden_data[start_index:end_index].strip()

        # reading atoms info
        atom_lines = fragment.split("\n")
        for line in atom_lines:
            split_line = line.split()
            if not len(split_line) == 6:
                raise Exception("Class Molecule. Unrecognized atom in molden file")

            position = [float(a)*0.529177249 for a in split_line[3:]]  # make a list of coordinates in Angstroms
            self.atoms.append(Atom(position=position, label=split_line[0], number=split_line[1]))

        self.update_molecule_limits()

    def update_molecule_limits(self):
        first_atom_podition = self.atoms[0].position
        
        self.x_min = self.x_max = first_atom_podition[0]
        self.y_min = self.y_max = first_atom_podition[1]
        self.z_min = self.z_max = first_atom_podition[2]

        for atom in self.atoms:
            self.x_min = min(atom.position[0], self.x_min)
            self.x_max = max(atom.position[0], self.x_max)
            self.y_min = min(atom.position[1], self.y_min)
            self.y_max = max(atom.position[1], self.y_max)
            self.z_min = min(atom.position[2], self.z_min)
            self.z_max = max(atom.position[2], self.z_max)


    def read_basis_functions_from_molden(self):
        # finding the [GTO] block in the file
        start_block = "[GTO]"
        end_block = "[5D]"
        start_index = self.molden_data.find(start_block) + len(start_block)
        if start_index == -1:
            raise Exception("Class Molecule. There is no [GTO] block in the file")

        end_index = self.molden_data.find(end_block, start_index)
        if end_index == -1:
            fragment = self.molden_data[start_index:]
        else:
            fragment = self.molden_data[start_index:end_index].strip()

        orbital_shells = ("s", "p", "d", "f")

        # reading gto info
        blocks = fragment.split("\n\n")  # dividing the GTO part into blocks coresponding to atoms
        # if not len(blocks) == len(self.atoms):
        #     raise Exception("Class Molecule. The number of GTO blocks is not equal to the number of atoms")
        for block in blocks:
            for GTO in re.split(r'(?=[spdf])', block)[1:]:  # dividing every block into fragments corresponding to GTOs
                GTO_lines = GTO.split("\n")
                shell = GTO_lines[0].split()[0]
                exponents = [float(a.split()[0]) for a in GTO_lines[1:] if a]
                coefficients = [float(a.split()[1]) for a in GTO_lines[1:] if a]
                number_of_atom = blocks.index(block)
                position = self.atoms[number_of_atom].position

                # print(shell, exponents, coefficients)

                if shell == "s":
                    self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position))
                elif shell == "p":
                    self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="x"))
                    self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="y"))
                    self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="z"))
                elif shell == "d":
                    self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="z2"))
                    self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="xz"))
                    self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="yz"))
                    self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="x2y2"))
                    self.basis_functions.append(BasisFunction(coefficients, exponents, shell, position, index="xy"))


            # i = 1
            # while i < len(lines):
            #     label = lines[i].split()[0]
            #     gto_numbers = lines[i].split()[1]
            #     for n in range(0, gto_numbers):
            #         exponent, coeff = lines

    def read_orbitals_from_molden(self):
        start_block = "[MO]"
        start_index = self.molden_data.find(start_block) + len(start_block)
        if start_index == -1:
            raise Exception("Class Molecule. There is no [MO] block in the file")

        fragment = self.molden_data[start_index:]
        blocks = re.split(r'(?=Sym=)', fragment)[1:]
        for block in blocks:
            MO_lines = [l.strip() for l in block.split('\n')]
            # print(MO_lines)
            occ = float(MO_lines[3].split("=")[1])
            # print(occ)
            coefs = [float(a.split()[1]) for a in MO_lines[4:] if a]
            # print(coefs)
            self.orbitals.append(Orbital(occ, coefs))


    def draw_molecule(self, ax):
        for atom in self.atoms:
            ax.scatter(atom.position[0], atom.position[1], atom.position[2], color=atom.ball_color, s=atom.ball_size, alpha=1)

    def scfp_density_at_point(self, point):
        if not self.density_matrix.size == 0:
            basis_functions_values = np.array([bf.value_at_point(point) for bf in self.basis_functions])

            # плотность ρ = φᵀ P φ
            return float(basis_functions_values @ self.density_matrix @ basis_functions_values)
        elif not len(self.orbitals) == 0:
            res = 0
            for orb in self.orbitals:
                if orb.occupancy == 0:
                    continue
                res += orb.occupancy*(self.orbital_value_at_point(point, self.orbitals.index(orb))**2)
            return res

    def orbital_value_at_point(self, point, number):
        basis_functions_values = np.array([bf.value_at_point(point) for bf in self.basis_functions])
        return(np.sum(basis_functions_values @ self.orbitals[number].coefficients))