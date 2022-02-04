from ase.collections import g2
from ase.build import molecule
import ase.io.vasp
from ase.visualize import view

"""
Now, using the molecule object. If you don't know a molecule is or not in a collection, use the line 11, which prints 
the name of any molecule in a collection.
"""

print(g2.names)  # Will print the name of all molecules in a collection.
atoms = molecule('PH3', vacuum=10)  # Create PH3 in a 10x10x10 box of vacuum.
view(atoms)

# Write a file called POSCAR containing features of our PH3 molecules.
ase.io.vasp.write_vasp('POSCAR', atoms, sort=True, direct=True)
