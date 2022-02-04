from ase import Atoms
from ase.visualize import view

"""
Start by using the simple atom object. The following script creates a PdHCl complex.
"""

atoms = Atoms('PdHCl', positions=[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [-2.1, 0.0, 0.0]])
atoms.center(vacuum=10)
view(atoms)
