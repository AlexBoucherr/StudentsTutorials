from ase.build import bulk, fcc100, fcc110, fcc111, fcc211, fcc111_root, surface
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.io import read
import ase.io.vasp

canh = read('/home/alexbou/CONTCAR_CaNH_isif')
slab = surface(canh, (0, 0, 1), layers=2, vacuum=5, periodic=True) * (2, 2, 1)

n_atoms = []
for atom in slab:
    if atom.symbol == 'N' and round(atom.position[2], 3) == 15.425:
        n_atoms.append(atom.index)
view(slab)

for atom in slab:
    if atom.index in n_atoms:
        atom.position[2] = 4.660
        atom.position[0] += 1.296
view(slab)

min_heigh = 100
for atom in slab:
    if atom.position[2] < min_heigh:
        min_heigh = atom.position[2]
view(slab)

for atom in slab:
    atom.position[2] = atom.position[2] - min_heigh

frozen_atoms = []
for atom in slab:
    if round(atom.position[2], 3) <= 3.464:
        frozen_atoms.append(atom.index)

slab.set_constraint(FixAtoms(indices=frozen_atoms))
view(slab)

ase.io.vasp.write_vasp('/home/alexbou/POSCAR_CaNH001_v2', slab, vasp5=True, sort=True, direct=True)