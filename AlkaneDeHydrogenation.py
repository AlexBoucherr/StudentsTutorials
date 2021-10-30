# First of all, we need to import all the 'tools' we will need to write ur script.
import ase
from ase import Atoms
from ase.build import surface, bulk, add_vacuum, molecule, add_adsorbate
from ase.visualize import view
from ase.constraints import FixAtoms
from ase.collections import g2
import numpy as np
import ase.io.vasp

# First of all, generate a simple bulk.
bulkMetal = bulk('Pd',  # The metal you're interested in.
                 a=3.929,  # Its lattice constant.
                 cubic=False)  # Use cubic=False for primitice cell.

# Want a bigger system? You can multiply any existant system by any 3D vector to make it bigger. For instance:
bulk2 = bulkMetal * (3, 3, 3)  # Generates a 27x times bigger bulk unit cell.

# Now, let's try and generate a slab out of our bulk material. To do so, use the function 'surface'.
slab = surface(bulkMetal,  # The initial system we use to generate our slab (our bulk).
               indices=(1, 1, 1),  # Miller indices.
               layers=5,  # The number of layers we want in the slab.
               vacuum=5,  # Put 5A vacuum on each side of the slab.
               periodic=True)  # Makes the system 3D periodic.

# Once again, if we want a larger system, multiply by a vector. For instance:
slab = slab * (5, 5, 1)
# You can use the following command to control the vacuum layer on top of your slab:
add_vacuum(slab, vacuum=5)  # Add 5 angstrom of vacuum on each side of the slab.
add_vacuum(slab, vacuum=-5)  # Remove 5 angstrom of vacuum on each side of the slab.

'''
So far we generated a 5x5 slab, 5 layers in Thickness. So as to make calculations faster, we can freeze multiple layers.
To do so, we can create what's called a function in Python, which is basically a procedure working as a mathematical
function, taking one or several argument and following a defined procedure.
'''


# Let's create a function that freezes a chosen number of layers in our slab.


def Freeze(slab, frozenlayers):
    # Create a list containing z-position of every layers.
    zCoordinates = []
    for atom in slab:
        if atom.position[2] not in zCoordinates:
            zCoordinates.append(atom.position[2])

    # Now, remove the layers on top. These layers will be relaxed upon relaxation of our system.
    if frozenlayers > len(zCoordinates):
        print('Error: Too many layers selected! Choose smaller layer parameter.')
        exit()  # If error, break the program.

    zCoordinates.sort()
    nRelax = len(zCoordinates) - frozenlayers  # Number of layers we will freeze.
    for i in range(0, nRelax):
        zCoordinates.pop(0)

    # Create a list containing the index of atoms we want to freeze.
    frozenAtoms = []
    for atom in slab:
        if atom.position[2] not in zCoordinates:
            frozenAtoms.append(atom.index)

    # Freeze them!
    slab.set_constraint(FixAtoms(frozenAtoms))

    # The line 'return' is used to provide the output of the function.
    return slab


# Now we simply call our function and specify the parameters to freeze the appropriate number of layers!
slab = Freeze(slab=slab, frozenlayers=2)

"""
We just figured out how to build bulk and slab. Now, let's figure out how to add an adsorbate on top of it.
The first step is to generate the adsorbate!
"""

# First of all, you can use 'Atoms' and set atomic positions yourself:
adsorbate = Atoms('C2H6',
                  positions=[[-1.45, 0.00, 0.00],  # C1, cartesian coordinates.
                             [1.45, 0.00, 0.00],  # C2
                             [2.22, 1.68, -1.00],  # H1
                             [2.22, 0.03, 1.96],  # H2
                             [2.22, -1.71, -0.95],  # H3
                             [-2.22, 1.71, 0.95],  # H4
                             [-2.22, -0.03, -1.96],  # H5
                             [-2.22, -1.68, 1.00]])  # H6

# view(molecule)
# There is (fortunatelly) a more convenient way to generate comon molecules with ASE!
adsorbate = molecule('C2H6')  # Convenient way to generate simple alkane.
adsorbate.center(vacuum=0)
view(adsorbate)


# What if you need to rotate the adsorbate? Create a 'rotation' function!
def Rotation(adsorbate, angle):
    adsorbate.center(vacuum=10)
    for atom in adsorbate:
        rotX = (atom.position[0] * np.sin(angle)) - (atom.position[1] * np.cos(angle))
        rotY = (atom.position[0] * np.sin(angle)) + (atom.position[1] * np.cos(angle))
        rotZ = atom.position[2]
        rotPosition = np.array([rotX, rotY, rotZ])
        atom.position = rotPosition
    adsorbate.center(vacuum=0)
    return adsorbate


adsorbate = Rotation(adsorbate=adsorbate, angle=1.0)
view(adsorbate)

# Generate POSCAR file.
ase.io.vasp.write_vasp('POSCAR_molecule',  # Name of output.
                       adsorbate,  # molecule.
                       direct=True,  # Keep it True!
                       vasp5=True,  # Keep it True!
                       sort=True, # Keep it True!
                       ignore_constraints=False)  # Keep it False!
