from ase.build import fcc100
from ase.constraints import FixAtoms, FixedLine
from ase.visualize import view

"""
The following large block of intructions is the 'Freezing' function. You don't need to be able to write such a function
yourself, you can simply copy-paste it in your own codes!
"""


def Freezing(file, relax_layers, partFreeze):
    if not partFreeze:
        # Create a list containing z_position of each atomic layers.
        z_positions = []
        for atom in file:
            if atom.position[2] not in z_positions:
                z_positions.append(atom.position[2])

        z_positions.sort()  # Sort the list.

        # Remove layers z_position depending on the number of slabs we need to freeze.
        if relax_layers != 0:
            # Remove the appropriate number of z_positions.
            for x in range(0, relax_layers):
                z_positions.remove(max(z_positions))
            file.set_constraint(FixAtoms(indices=[i.index for i in file if i.position[2] in z_positions]))
        return file

    if partFreeze:
        # Create a list containing z_position of each atomic layers.
        z_positions = []
        for atom in file:
            if atom.position[2] not in z_positions:
                z_positions.append(atom.position[2])

        z_positions.sort()  # Sort the list.

        # Remove layers z_position depending on the number of slabs we need to freeze.
        if relax_layers != 0:
            # Remove the appropriate number of z_positions.
            for x in range(0, relax_layers):
                z_positions.remove(max(z_positions))

        indicesTotFeeze = []
        indicePartFreeze = []
        for atom in file:
            if atom.position[2] in z_positions:
                indicesTotFeeze.append(atom.index)
            if atom.position[2] not in z_positions:
                indicePartFreeze.append(atom.index)

        constraint1 = FixAtoms(indices=indicesTotFeeze)
        constraint2 = FixedLine(indicePartFreeze, direction=[0, 0, 1])
        file.set_constraint([constraint1, constraint2])
        return file


slab = fcc100('Pd', size=(3, 3, 5), periodic=True, vacuum=10)
slab = Freezing(file=slab, relax_layers=2, partFreeze=False)
view(slab)