import numpy as np

class Orbital:
    def __init__(self, occupancy, coefficients):
        self.occupancy = occupancy
        self.coefficients = np.array(coefficients)

    def __repr__(self):
        return (f"""Occupancy = {self.occupancy}
                coefficients = {self.coefficients}""")