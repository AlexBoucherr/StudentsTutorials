from ase.build import bulk, fcc100, fcc110, fcc111, fcc211, fcc111_root, surface
from ase.visualize import view

"""
Create a slab using some of the ASE functions: fcc111, fcc100, fcc110, hcp0001... Or the surface function.
"""
# Option 1: a 100 slab of Pd.
slab = fcc100('Pd', size=(3, 3, 5), periodic=True, vacuum=10)
view(slab)

# Option2: Use the surface function.
material = bulk('Pd', crystalstructure='fcc', cubic=False)
slab = surface(material, (6, 1, 1), layers=20, vacuum=10, periodic=True) * (1, 5, 1)
view(slab)