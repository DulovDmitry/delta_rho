class Atom:
	def __init__(self, position, label="", number=0, charge=0):
		self.position = position
		self.label = label
		self.number = number
		self.charge = charge

		colors_dict = {
		"H" : "red", #"#cccccc",
		"C" : "#666666",
		"N" : "#0066ff",
		}

		sizes_dict = {
		"H" : 25,
		"C" : 35,
		"N" : 35,
		}

		self.ball_color = colors_dict[label]
		self.ball_size = sizes_dict[label]

	def __repr__(self):
		return (f"""Atom label = {self.label}
				number = {self.number}
				position = {self.position}""")
