from ase.io import read  # Use to read files.
from ase.constraints import FixAtoms  # Use to freeze atoms.
from ase.visualize import view  # Use to view objects
from ase import Atoms  # Use to createatoms objects
import ase.io.vasp  # Use to write files.
from ase.build import add_adsorbate  # Add an adsorbate on a surface

# Calls the read function to read the desired file. Set the path (between brackets) to desired file.
silica = read('/home/alex/Downloads/BareSilica/Pristine_fullH_001/Tier2/CONTCAR') * (1, 1, 1)  # Controls size of cells.

# Remove existing constraints.
silica.constraints = False
frozen_atoms = []

# Run over each atom in the system. If they are below a given height (13.0 Angstrom) we will freeze them.
for atom in silica:
    if atom.position[2] < 13.000:
        frozen_atoms.append(atom.index)

# Apply a constraint freezing desired atoms.
silica.set_constraint(FixAtoms(indices=frozen_atoms))

# Read the complex. Use sdf files from PubChem.
w_complex = read('/home/alex/Downloads/Structure2D_CID_16217328.sdf')

# Put the complex in a 10x10x10 Angstrom box of vacuum.
w_complex.center(vacuum=10)
view(w_complex)

# Creating adsorption group.
del_atoms = [40, 39, 38, 37, 36, 35, 14, 13, 1]  # List of atoms part of the leaving group.
del_atoms.sort(reverse=True)
adsorb = True

# If adsorb = True, remove the leaving group. Use this structure to put on the silica surface.
if adsorb:
    for x in del_atoms:
        del w_complex[x]
    view(w_complex)

departing = False
# If departing = True, only the atoms of the departing group are left. DONT SET DEPARTING AND ADSORB TRUE BOTH AT ONCE!
if departing:
    del_atoms2 = [atom.index for atom in w_complex if atom.index not in del_atoms]
    del_atoms2.sort(reverse=True)
    for x in del_atoms2:
        del w_complex[x]
    # Create a single atom of Hydrogen and adsorb it to the leaving group on N to balance charges.
    h = Atoms('H')
    add_adsorbate(w_complex, h, 0, position=[13.1, 12])
    view(w_complex)

# Add adsorbate on top of silica. Start by deleting a H atom.
del silica[6]

# Add the W-group obtained by setting adsrob = True.
add_adsorbate(silica, w_complex, 2.5, position=[5.5, 3.7])
view(silica)

# Set the name of the file created as desired (first argument).
ase.io.vasp.write_vasp('POSCAR_name', w_complex, direct=True, vasp5=True)
