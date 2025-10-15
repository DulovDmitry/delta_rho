import numpy as np

ANGSTROM_TO_BOHR = 1.8897259886
ANGSTROM_TO_BOHR_SQUARED = 1.8897259886**2

class BasisFunction:
    def __init__(self, coefficients, exponents, shell, position, label=None, index=None):
        self.coefficients = np.array(coefficients)
        self.exponents = np.array(exponents)
        self.shell = shell
        self.position = position
        self.label = label
        self.index = index

        self.s_res_dict = dict()

    def __repr__(self):
        return (f"""Function label = {self.label}
                position = {self.position}
                coefficients = {self.coefficients}
                exponents = {self.exponents}
                shell = {self.shell}
                index = {self.index}""")

    def value_at_point(self, point):
        if self.shell == "s":
            return self.s(point)
        elif self.shell == "p":
            return self.p(point)
        elif self.shell == "d":
            return self.d(point)
        else:
            raise Exception("Class BasisFunction. Unknown orbital shell")

    def s(self, point):
        squared_radius_vector = ANGSTROM_TO_BOHR_SQUARED*((self.position[0] - point[0])**2 + (self.position[1] - point[1])**2 + (self.position[2] - point[2])**2)
        if squared_radius_vector > 10:
            return 0

        norm_coefs = np.pow(2*self.exponents/np.pi, 3/4)
        # res1 = norm_coefs*self.coefficients*np.exp(-self.exponents*squared_radius_vector)
        if squared_radius_vector in self.s_res_dict.keys():
            res = self.s_res_dict[squared_radius_vector]
        else:
            res = np.sum(self.coefficients*np.exp(-self.exponents*squared_radius_vector))
            # res=0
            # for coef, exp in zip(self.coefficients,self.exponents):
            #     res += coef*np.exp(-exp*squared_radius_vector)
            self.s_res_dict[squared_radius_vector] = res

        return res

    def p(self, point):
        # return 0
        squared_radius_vector = ANGSTROM_TO_BOHR_SQUARED*((self.position[0] - point[0])**2 + (self.position[1] - point[1])**2 + (self.position[2] - point[2])**2)
        if squared_radius_vector > 10:
            return 0
        # norm_coefs = np.pow(128*(self.exponents**5)/(np.pi**3), 1/4)
        if self.index == "x":
            res = self.coefficients*ANGSTROM_TO_BOHR*(point[0] - self.position[0])*np.exp(-self.exponents*squared_radius_vector)
        elif self.index == "y":
            res = self.coefficients*ANGSTROM_TO_BOHR*(point[1] - self.position[1])*np.exp(-self.exponents*squared_radius_vector)
        elif self.index == "z":
            res = self.coefficients*ANGSTROM_TO_BOHR*(point[2] - self.position[2])*np.exp(-self.exponents*squared_radius_vector)
        else:
            raise Exception("Class BasisFunction. Unknown p-orbital")
        return np.sum(res)

    def d(self, point):
        return 0
