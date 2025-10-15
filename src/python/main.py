import numpy as np
import matplotlib.pyplot as plt

from molecule import Molecule
from grid import SphericalGrid, RegularOrthogonalGrid, IrregularOrthogonalGrid


json_file_name = "D:/lambda_project/dihydrogen_opt.json"
# json_file_name = "D:/lambda_project/BSF-33_opt_b3lyp.json"
# molden_file_name = "D:/lambda_project/dihydrogen_opt.molden.input"
molden_file_name = "D:/lambda_project/BSF-33_opt_b3lyp.molden.input"

mol2 = Molecule.from_molden(molden_file_name)
# print(mol2.orbitals)
# print(mol2.basis_functions)


# mol = Molecule.from_json(json_file_name)


center_point = (0,0,0)
# center_point = (mol.x_coordinates[0],mol.y_coordinates[0],mol.z_coordinates[0])

# grid = SphericalGrid(center = mol2.atoms[0].position)
# grid = RegularOrthogonalGrid(center = mol2.atoms[0].position, x_points=50, y_points=50, z_points=25)
grid = RegularOrthogonalGrid.from_molecule(mol2,20,20,20)
# grid = IrregularOrthogonalGrid(center = mol2.atoms[0].position)
values = grid.calculate_scfp_values(mol2)
print(grid.scfp_integral())
print(mol2.scfp_density_at_point(mol2.atoms[0].position))


print(f"Shape of x, y, z: {grid.x.shape}")
print(f"Min and max values - x: ({grid.x.min():.2f}, {grid.x.max():.2f}), y: ({grid.y.min():.2f}, {grid.y.max():.2f}), z: ({grid.z.min():.2f}, {grid.z.max():.2f})")
print(f"Center point: ", grid.center_point)
print(f"Density: max value = {max(values)}, min value = {min(values)}")

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
# grid.draw_isosurface(ax, 0.25)s
# grid.draw_grid_points(ax, molecule=mol2, type="scfp")
# grid.draw_grid_points(ax)
grid.draw_grid_points(ax, molecule=mol2, type="orbital", number=100)
mol2.draw_molecule(ax)
ax.set_box_aspect((1,1,1))
plt.show()