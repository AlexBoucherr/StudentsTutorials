########################################################################################################################
# ASE learning script 2. → Tutorial to create vaccancies & adatoms on perfect surfaces.                                #
# Version 1.0 created on 21/02/2021 by Alex.                                                                           #
# Last update: 01/06/2021.                                                                                             #
########################################################################################################################
"""
Last time, we learnt how to create various slabs using two tools, 'Bulk' and 'Surface'. Today, we are using a new tool,
'fcc111', allowing you to directly create a (111) surface without creating the bulk at first. On this (111) surface,
we will learn how to create vacancies and add-atoms (two common defects in materials' surfaces).
"""
# First of all, we import all the tools we will need:
from ase.build import fcc111, fcc100, fcc110, add_adsorbate
from ase import Atoms
from ase.constraints import FixAtoms, FixedLine
import ase.io.vasp
from ase.visualize import view
from ase.io import read
import numpy as np

# Then, we create all the global variables we will need:
READ_FILE = False                                     # If needed to read external file. Read+path or False.
LAYERS = 5                                            # Select the number of layers you want.
SURFACE = 111                                         # Nature of the surface to generate.
VAC = 15/2                                            # Select thickness of Vacuum layer.
SIZE = (5, 5, LAYERS)                                 # The size of the system.
NAME = 'POSCAR_100_1Ad'                               # Name the file (cf. previous tutorials).
LAT = 3.929                                           # Lattice constant (Pd), RPBE-D3 ISIF calculations.
RELAXED = 2                                           # Relaxed layers (before modification) (0, 1, 2 or 'all').

DEL = 0                                               # Delete an Atom? 0, 1, 2 or 3.
ADD_AT = 0                                            # Numbers of atoms to add (0, 1, 2 or 3).
GEOMETRY = 'Pyr'                                      # Pyr, Sq, Sq_shift.
low_CN = False                                        # False, 'CN1' or 'CN2'.

########################################################################################################################
# A few functions we will need.                                                                                        #
########################################################################################################################


def Freezing(slab, relaxed_layers):
    Z_coord = []
    for atom in slab:
        if atom.position[2] not in Z_coord:
            Z_coord.append(atom.position[2])

    Z_coord.sort()
    for x in range(0, relaxed_layers):
        Z_coord.remove(max(Z_coord))

    constraint = FixAtoms(indices=[i.index for i in PdSurf if i.position[2] in Z_coord])

    return constraint


'''
We now use a new tool, 'fcc111' to create the slab:
'''
if SURFACE == 111 and not READ_FILE:
    PdSurf = fcc111('Pd',                                 # The nature of the atom.
                    size=SIZE,                            # The size of the system.
                    a=LAT,                                # The lattice constant of Pd.
                    vacuum=VAC)                           # Bilayer of vacuum (on top & bottom as previously seen).

if SURFACE == 100 and not READ_FILE:
    PdSurf = fcc100('Pd',                                 # The nature of the atom.
                    size=SIZE,                            # The size of the system.
                    a=LAT,                                # The lattice constant of Pd.
                    vacuum=VAC)                           # Bilayer of vacuum (on top & bottom as previously seen).

if SURFACE == 110 and not READ_FILE:
    PdSurf = fcc110('Pd',                                 # The nature of the atom.
                    size=SIZE,                            # The size of the system.
                    a=LAT,                                # The lattice constant of Pd.
                    vacuum=VAC)                           # Bilayer of vacuum (on top & bottom as previously seen).

if READ_FILE:
    PdSurf = READ_FILE

# This command is more easy to use than the 'surface' tool introduced in the previous session.
# However, you can only create (111) slabs using the 'fcc111' tool.
'''
We will now put the slab on the bottom of the unit cell.
The procedure is exactly the same as the one used in the previous training.
'''
# This loop shift the z-component of the position of each atom.
if not READ_FILE:
    PdSurf.set_positions([[i.position[0], i.position[1], i.position[2]-VAC] for i in PdSurf])

'''
Now, we will freeze the inner layers. Once again the procedure is the same as used in previous training.
We 'FixAtoms' any atom on the inner layer. They are deeper enough so that the presence of the vacuum layer on top has no
effect on them.
'''
# If relaxed = all, we don't do anything.
if RELAXED == 'all':
    pass

# Otherwise, use the freezing function to generate a given number of frozen layers.
if RELAXED != 'all':
    PdSurf.set_constraint(Freezing(PdSurf, RELAXED))

########################################################################################################################
# BELOW: INSTRUCTIONS TO CREATE UP TO 2 VACANCIES.                                                                     #
########################################################################################################################
'''
To delete atoms in a structure, we use the generic command 'del'. For your project, you will have to create single or
double vacancies. Let's create a algorithm, using loops and conditions, (these structures are fundamental in Python)
to define which atom(s) we want to delete.
'''

# If yes, then show him the index of the atoms on the top layer
if DEL != 0:
    view(PdSurf)
    Z_coord2 = []
    for k in PdSurf:
        if k.position[2] not in Z_coord2:
            Z_coord2.append(k.position[2])
    z_max = max(Z_coord2)
    for i in PdSurf:
        if i.position[2] == z_max:
            print(i.index, [i.position[0], i.position[1]])
        else:
            pass

# Now, two situations. Delete 1 atom or delete 2 atoms?
if DEL == 1:                                                                    # If previous answer = 1.
    DelPd1 = int(input('Choose the index of the Pd atom to delete:\n'))           # Ask index of the atom to delete.
    del PdSurf[[atom.index for atom in PdSurf if atom.index == DelPd1]]           # The atom is deleted.
elif DEL == 2:
    L = []                                                                        # Create a list to store indices.
    DelPd2_1 = int(input('Choose the index of the first Pd atom to delete:\n'))   # Ask for first index.
    L.append(DelPd2_1)                                                            # Put it in the list.
    DelPd2_2 = int(input('Choose the index of the second Pd atom to delete:\n'))  # Ask for the second index.
    L.append(DelPd2_2)                                                            # Put it in the list.
    L.sort()                                                                      # Sort the list.
    del PdSurf[[atom.index for atom in PdSurf if atom.index == L[1]]]             # We delete big_index first.
    del PdSurf[[atom.index for atom in PdSurf if atom.index == L[0]]]             # And then small_index.
elif DEL == 3:
    L = []
    DelPd2_1 = int(input('Choose the index of the first Pd atom to delete:\n'))
    L.append(DelPd2_1)
    DelPd2_2 = int(input('Choose the index of the second Pd atom to delete:\n'))
    L.append(DelPd2_2)
    DelPd2_3 = int(input('Choose the index of the third Pd atom to delete:\n'))
    L.append(DelPd2_3)
    L.sort()
    del PdSurf[[atom.index for atom in PdSurf if atom.index == L[2]]]
    del PdSurf[[atom.index for atom in PdSurf if atom.index == L[1]]]
    del PdSurf[[atom.index for atom in PdSurf if atom.index == L[0]]]

if DEL == 0:
    pass

########################################################################################################################
# BELOW: INSTRUCTIONS TO GENERATE ADD-ATOMS.                                                                           #
########################################################################################################################
'''
first things first, let's create the objects you will insert on the Pd surface. 
For your project you will insert up to 3 ad-atoms.Let's first create those structures.
'''
Pd_Pd = 2.772                                           # Pd-Pd distance.
h = 2.1                                                 # The height to insert the ad-atoms.

Pd1 = Atoms('Pd')                                       # We create a single Pd atom.
Pd2 = Atoms('2Pd',                                      # We create 2 Pd atoms.
            positions=[[0, 0, 0],                       # We specify the positions of the 2 atoms.
                       [Pd_Pd, 0, 0]])
Pd3 = Atoms('3Pd',                                      # We create 3 Pd atoms.
            positions=[[0, 0, 0],                       # We specify the position of each atom.
                       [Pd_Pd, 0, 0],
                       [2*Pd_Pd, 0, 0]])
Pd4 = Atoms('4Pd',
            positions=[[0, 0, 0],
                       [Pd_Pd, 0, 0],
                       [0, Pd_Pd, 0],
                       [Pd_Pd, Pd_Pd, 0]])

Pd4_pyr = Atoms('4Pd',
                positions=[[0, 0, 0],
                           [Pd_Pd, 0, 0],
                           [Pd_Pd / 2, 0.7 * Pd_Pd, 0],
                           [Pd_Pd / 2, 0.35 * Pd_Pd, 0.7 * Pd_Pd]])

Pd4_Sq2 = Atoms('4Pd',
                positions=[[0, 0, 0],
                           [Pd_Pd, 0, 0],
                           [0.5 * Pd_Pd, Pd_Pd, 0],
                           [1.5 * Pd_Pd, Pd_Pd, 0]])

Pd5 = Atoms('5Pd',
            positions=[[0, 0, 0],
                       [Pd_Pd, 0, 0],
                       [0, Pd_Pd, 0],
                       [Pd_Pd, Pd_Pd, 0],
                       [0.5 * Pd_Pd, 0.5 * Pd_Pd, 0.7 * Pd_Pd]])

# We create different conditions depending on the previous answer.
# To insert atoms on our slabs, we use the ASE.build tool 'add_adsorbate':
if ADD_AT == 1:                                    # If answer = 1. Adsorbate = Single Pd atom.
    if not low_CN:
        AdAtoms = Pd1
        X = float((PdSurf.cell[0][0] / 2) + 0.77 * (PdSurf.cell[1][0] / 2))  # Adsorbate x-position.
        Y = float(PdSurf.cell[1][1] / 2)                                     # Adsorbate y-position.
        add_adsorbate(PdSurf,                                                # The Slab.
                      AdAtoms,                                               # The adsorbate.
                      h,                                                     # z-position above the slab.
                      position=[X, Y])                                       # Adsorbate x & y-positions.
    elif low_CN == 'CN2' and SURFACE == 100:
        AdAtoms = Pd1
        X = float((PdSurf.cell[0][0] / 2) + 0.77 * (PdSurf.cell[1][0] / 2)) - (2.778 / 2)
        Y = float(PdSurf.cell[1][1] / 2)
        add_adsorbate(PdSurf,
                      AdAtoms,
                      h,
                      position=[X, Y])
    elif low_CN == 'CN1' and SURFACE == 100:
        AdAtoms = Pd1
        X = float((PdSurf.cell[0][0] / 2) + 0.77 * (PdSurf.cell[1][0] / 2)) - (2.778 / 2)
        Y = float(PdSurf.cell[1][1] / 2) - (2.778 / 2)
        add_adsorbate(PdSurf,
                      AdAtoms,
                      h,
                      position=[X, Y])
    elif low_CN == 'CN1' and SURFACE == 111:
        view(PdSurf)
        z_max = 0
        for atom in PdSurf:
            if atom.position[2] > z_max:
                z_max = atom.position[2]

        site = int(input('Select site: \n'))
        X = [atom.position[0] for atom in PdSurf if atom.index == site]
        X = X[0]
        Y = [atom.position[1] for atom in PdSurf if atom.index == site]
        Y = Y[0]

        AdAtoms = Pd1
        add_adsorbate(PdSurf,
                      AdAtoms,
                      h + 0.3,
                      position=[X, Y])

        del PdSurf.constraints  # Delete existing constraints.

        constraint1 = Freezing(PdSurf, RELAXED + 1)  # +1 for the extra layer accounted from the adatom!
        constraint2 = FixedLine(a=[atom.index for atom in PdSurf if atom.position[2] > z_max],
                                direction=[0, 0, 1])

        PdSurf.set_constraint([constraint1, constraint2])

        view(PdSurf)
    elif low_CN == 'CN2' and SURFACE == 111:
        view(PdSurf)
        z_max = 0
        for atom in PdSurf:
            if atom.position[2] > z_max:
                z_max = atom.position[2]

        site1 = int(input('Select site 1: \n'))
        site2 = int(input('select site 2 (adjacent to site 1): \n'))

        for atom in PdSurf:
            if atom.index == site1:
                site1 = atom
            if atom.index == site2:
                site2 = atom
        X = site1.position[0] + ((site2.position[0] - site1.position[0]) / 2)
        Y = site1.position[1] + ((site2.position[1] - site1.position[1]) / 2)

        AdAtoms = Pd1
        add_adsorbate(PdSurf,
                      AdAtoms,
                      h,
                      position=[X, Y])

        del PdSurf.constraints  # Delete existing constraints.

        constraint1 = Freezing(PdSurf, RELAXED + 1)  # +1 for the extra layer accounted from the adatom!
        constraint2 = FixedLine(a=[atom.index for atom in PdSurf if atom.position[2] > z_max],
                                direction=[0, 0, 1])

        PdSurf.set_constraint([constraint1, constraint2])
        view(PdSurf)

# The position of the atom on the slab is obsviously not optimal. However, it should a good enough guess for
# a geometry optimisation process!

if ADD_AT == 2:
    AdAtoms = Pd2
    X = float((PdSurf.cell[0][0] / 2) + 0.38 * (PdSurf.cell[1][0] / 2))  # We will use it to position the adsorbate.
    Y = float(PdSurf.cell[1][1] / 2)  # We will use it to position the adsorbate.
    add_adsorbate(PdSurf,
                  AdAtoms,
                  h,
                  position=[X, Y])

if ADD_AT == 3 and SURFACE == 111:
    AdAtoms = Pd3
    X = float((PdSurf.cell[0][0] / 2) + 0.02 * (PdSurf.cell[1][0] / 2))  # We will use it to position the adsorbate.
    Y = float(PdSurf.cell[1][1] / 2)  # We will use it to position the adsorbate.
    add_adsorbate(PdSurf,
                  AdAtoms,
                  h,
                  position=[X, Y])

if ADD_AT == 3 and SURFACE == 100:
    AdAtoms = Pd3
    X = float((PdSurf.cell[0][0] / 2) + 0.02 * (PdSurf.cell[1][0] / 2))  # We will use it to position the adsorbate.
    Y = float(PdSurf.cell[1][1] / 2)  # We will use it to position the adsorbate.
    add_adsorbate(PdSurf,
                  AdAtoms,
                  h,
                  position=[X - 2.778, Y])

if ADD_AT == 4:
    if GEOMETRY == 'Sq':
        AdAtoms = Pd4
        X = float((PdSurf.cell[0][0] / 2) + 0.38 * (PdSurf.cell[1][0] / 2)) - Pd_Pd  # We will use it to position the adsorbate.
        Y = float(PdSurf.cell[1][1] / 2) - Pd_Pd  # We will use it to position the adsorbate.
        add_adsorbate(PdSurf,
                      AdAtoms,
                      h,
                      position=[X, Y])
    if GEOMETRY == 'Pyr':
        AdAtoms = Pd4_pyr
        X = float((PdSurf.cell[0][0] / 2) + 0.38 * (
                    PdSurf.cell[1][0] / 2)) - Pd_Pd  # We will use it to position the adsorbate.
        Y = float(PdSurf.cell[1][1] / 2) - Pd_Pd  # We will use it to position the adsorbate.
        add_adsorbate(PdSurf,
                      AdAtoms,
                      h,
                      position=[X, Y])
    if GEOMETRY == 'Sq_shift':
        AdAtoms = Pd4_Sq2
        X = float((PdSurf.cell[0][0] / 2) + 0.38 * (
                PdSurf.cell[1][0] / 2)) - Pd_Pd  # We will use it to position the adsorbate.
        Y = float(PdSurf.cell[1][1] / 2) - Pd_Pd  # We will use it to position the adsorbate.
        add_adsorbate(PdSurf,
                      AdAtoms,
                      h,
                      position=[X, Y])

if ADD_AT == 5:
    AdAtoms = Pd5
    X = float((PdSurf.cell[0][0] / 2) + 0.38 * (PdSurf.cell[1][0] / 2)) - Pd_Pd  # We will use it to position the adsorbate.
    Y = float(PdSurf.cell[1][1] / 2) - Pd_Pd  # We will use it to position the adsorbate.
    add_adsorbate(PdSurf,
                  AdAtoms,
                  h,
                  position=[X, Y])

if ADD_AT == 0:                                    # If previous answer is 0, we don't do anything.
    pass


view(PdSurf)
'''
Now, all we need is to create a VASP-compatible file with the structure we just created.
The procedure is the same as the one previously used in the first tutorial.
'''
# We create a VASP-compatible file:
ase.io.vasp.write_vasp(NAME,                            # The name of the file we create.
                       PdSurf,                          # The object we use to create the file (our surface).
                       direct=True,
                       vasp5=True,                      # VASP5 compatible file.
                       sort=True,
                       ignore_constraints=False)        # If True, the 'FixAtoms set-up' will be ignored.

