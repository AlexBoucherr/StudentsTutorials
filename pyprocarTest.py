from ase import Atoms
import pyprocar
from ase.dft.kpoints import bandpath
from ase.visualize import view
import ase.io.vasp

atoms = Atoms('Pd2', positions=[[0, 0, 0], [2.7, 0, 0]])
atoms.center(vacuum=6)
ase.io.vasp.write_vasp('POSCAR', atoms, direct=True, vasp5=True)

pyprocar.kpath('POSCAR', 'KPOINTS', 10)
