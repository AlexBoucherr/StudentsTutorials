from ase.constraints import FixAtoms
from ase.io import read
from ase.calculators.vasp import Vasp2

slab = read('./CONTCAR')

# Remove existing constraints.
slab.constraints = False

# Create new constraints.
frozen_atoms = []
for atom in slab:
    if atom.symbol != 'Pd':
        frozen_atoms.append(atom.index)

slab.set_constraint(FixAtoms(indices=frozen_atoms))

# Create the calculator
dft = Vasp2(atoms=slab, directory='./', ncore=1, prec='normal', xc='rpbe', lwave=False, lcharg=False, ivdw=11,
            icharg=0, lmaxmix=2, ispin=2, kspacing=0.2, gamma=True, charge=0, encut=500, ediff=1E-9, lreal='auto',
            nelm=1500, algo='conjugate', lasph=True, laechg=False, lelf=False, ldipol=False, idipol=3, ismear=0,
            sigma=0.1, ediffg=1E-8, ibrion=1, nsw=500, potim=0.025)

# Attach the calculator to the atoms.
slab.calc = dft

# Get the energy forces and frequencies
slab.get_potential_energy()
