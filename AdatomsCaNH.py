from ase.io import read
from ase.visualize import view
from ase.build import molecule, add_adsorbate
import ase.io.vasp

slab = read('/home/alexbou/CONTCAR')
adsorbate = molecule('NH3')

"""
metal postion: [7.294, 4.996]
Ca top: [4.583, 2.498]
NH top: [4.281, 4.996]
"""
add_adsorbate(slab=slab, adsorbate=adsorbate, height=2.000, position=[4.281, 4.996])

ase.io.vasp.write_vasp('/home/alexbou/POSCAR_CaNHCo_NH3', slab, direct=True, sort=True, vasp5=True)
