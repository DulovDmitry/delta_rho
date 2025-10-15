import numpy as np
from molecule import Molecule
from scipy.integrate import simpson

class Grid:
	def __init__(self):
		self.scfp_values = list()

	def calculate_scfp_values(self, molecule: Molecule):
		if (not molecule.density_matrix.size == 0) or (not len(molecule.orbitals) == 0):
			for x, y, z in zip(self.x, self.y, self.z):
				self.scfp_values.append(molecule.scfp_density_at_point((x,y,z)))

			self.scfp_values = np.array(self.scfp_values)
			self.scfp_values_normalized = self.scfp_values/self.scfp_values.max()
			self.grid_points_sizes = self.scfp_values_normalized*20
		else:
			raise Exception("The molecule has niether density matrix nor orbitals")
		return self.scfp_values


	def draw_grid_points(self, ax, type=None, molecule=None, number=None):
		if type == "scfp":
			if len(self.scfp_values):
				ax.scatter(self.x, self.y, self.z, cmap="Blues", c=self.scfp_values, s=self.grid_points_sizes, alpha=0.4)
			else:
				raise Exception("You should first calculate scfp values using calculate_scfp_values() function")
		elif type == "orbital":
			if number == None:
				raise Exception("You should specify the number of orbital")

			values = list()
			for x, y, z in zip(self.x, self.y, self.z):
				values.append(molecule.orbital_value_at_point((x,y,z), number))

			ax.scatter(self.x, self.y, self.z, cmap="PiYG", c=values, s=7, alpha=0.4)
		else:
			ax.scatter(self.x, self.y, self.z, c="blue", s=5, alpha=0.4)

	def draw_isosurface(self, ax, value):
		x = np.asarray(self.x)
		y = np.asarray(self.y)
		z = np.asarray(self.z)

		tol = 0.01
		mask = self.scfp_values >= value

		ax.scatter(self.x[mask], self.y[mask], self.z[mask], c="red", s=15, alpha=0.5)

	def scfp_integral(self):
		raise NotImplementedError("Subclass must implement abstract method")


class SphericalGrid(Grid):
	def __init__(self, phi_points=20, theta_points=10, r_points=10, radius=3, center=(0,0,0)):
		super().__init__()
		self.delta_phi = 2*np.pi/phi_points
		phi_initial = self.delta_phi/2
		phi_final = 2*np.pi-(self.delta_phi/2)

		self.delta_theta = np.pi/theta_points
		theta_initial = self.delta_theta/2
		theta_final = np.pi-(self.delta_theta/2)

		self.delta_r = 2*radius/(2*r_points-1)
		r_initial = self.delta_r/2

		phi = np.linspace(phi_initial, phi_final, phi_points)
		theta = np.linspace(theta_initial, theta_final, theta_points)
		r = np.linspace(r_initial, radius, r_points)

		r, theta, phi = np.meshgrid(r, theta, phi) 

		# Convert to Cartesian coordinates
		x = center[0] + r * np.sin(theta) * np.cos(phi)
		y = center[1] + r * np.sin(theta) * np.sin(phi)
		z = center[2] + r * np.cos(theta)

		self.x = x.flatten()
		self.y = y.flatten()
		self.z = z.flatten()

		self.r = r.flatten()
		self.phi = phi.flatten()
		self.theta = theta.flatten()

		self.center_point = center

	
class RegularOrthogonalGrid(Grid):
	def __init__(self, x_length=4, y_length=4, z_length=4, x_points=10, y_points=10, z_points=10, center=(0,0,0)):
		super().__init__()
		self.x_points = np.linspace(center[0] - x_length/2, center[0] + x_length/2, x_points)
		self.y_points = np.linspace(center[1] - y_length/2, center[1] + y_length/2, y_points)
		self.z_points = np.linspace(center[2] - z_length/2, center[2] + z_length/2, z_points)

		x, y, z = np.meshgrid(self.x_points, self.y_points, self.z_points)

		self.x = x.flatten()
		self.y = y.flatten()
		self.z = z.flatten()

		self.delta_x = x_length/(x_points-1)
		self.delta_y = y_length/(y_points-1)
		self.delta_z = z_length/(z_points-1)
		self.delta_V = self.delta_x * self.delta_y * self.delta_z

		self.center_point = center

	@classmethod
	def from_molecule(cls, molecule: Molecule, x_points=20, y_points=20, z_points=20):
		x_margin = 2  # value in Angstroms
		y_margin = 2 
		z_margin = 2

		x_length = molecule.x_max - molecule.x_min + x_margin
		y_length = molecule.y_max - molecule.y_min + y_margin
		z_length = molecule.z_max - molecule.z_min + z_margin

		center_point = ((molecule.x_min+molecule.x_max)/2, (molecule.y_min+molecule.y_max)/2, (molecule.z_min+molecule.z_max)/2)

		return cls(x_length, y_length, z_length, x_points, y_points, z_points, center_point)

	def scfp_integral(self):
		# return np.sum(self.scfp_values) * self.delta_V
		Nx, Ny, Nz = len(self.x_points), len(self.y_points), len(self.z_points)
		f = self.scfp_values.reshape((Nx, Ny, Nz))
		I = simpson(
		        simpson(
		            simpson(f, self.z_points, axis=2),
		            self.y_points, axis=1),
		        self.x_points, axis=0
		    )
		return I


class IrregularOrthogonalGrid(Grid):
	def __init__(self, delta=0.01, scaling_factor=1.3, points=20, center=(0,0,0)):
		super().__init__()
		template = np.cumsum(((scaling_factor**np.linspace(0, points-1, points))*delta))
		print(template)

		self.x_points = np.append(center[0] + template, center[0] - template)
		self.y_points = np.append(center[1] + template, center[1] - template)
		self.z_points = np.append(center[2] + template, center[2] - template)

		self.x_points = np.append(center[0], self.x_points)
		self.y_points = np.append(center[1], self.y_points)
		self.z_points = np.append(center[2], self.z_points)
		
		print(self.x_points)

		self.x_points = np.sort(self.x_points)
		self.y_points = np.sort(self.y_points)
		self.z_points = np.sort(self.z_points)

		print(self.x_points)


		x, y, z = np.meshgrid(self.x_points, self.y_points, self.z_points)

		self.x = x.flatten()
		self.y = y.flatten()
		self.z = z.flatten()

		print(self.x)
		print(self.y)
		print(self.z)

		self.center_point = center

	def scfp_integral(self):
		# return np.sum(self.scfp_values) * self.delta_V
		Nx, Ny, Nz = len(self.x_points), len(self.y_points), len(self.z_points)
		f = self.scfp_values.reshape((Nx, Ny, Nz))
		I = simpson(
		        simpson(
		            simpson(f, self.z_points, axis=2),
		            self.y_points, axis=1),
		        self.x_points, axis=0
		    )
		return I
