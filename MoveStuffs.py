from ase.build import molecule, fcc111, add_adsorbate
from ase.visualize import view

slab = fcc111('Pd', size=(5, 5, 5), vacuum=10, periodic=True)
water = molecule('H2O')
for atom in water:
    atom.position *= -1
add_adsorbate(slab, water, position=[9.627, 7.146], height=2.5)

"""
vector [a1, a2] and [a3, a4], par=angles set given angle between the two vectors.
 indices = indices of the atoms to rotate.
"""
slab.set_dihedral(a1=122, a2=117, a3=125, a4=126, angle=125, indices=[125, 126, 127])
"""
Set the angle between atoms a1, a2 and a3. add=True, angle is added to existing angle. add=False: Sets a new angle.
indices = indices of atoms to rotate along with a3.
"""
slab.set_angle(a1=118, a2=117, a3=125, angle=20, add=True, indices=[125, 126, 127])
view(slab)