from ase.build import fcc100
from ase.visualize import view

"""
Remove an atom. Once again, you don't necessarily have to be able to build this function yourself!
"""


def Remove(file):
    view(file)
    indices = input('Enter the index of the atom(s) to delete:').split()

    atomsRemove = []

    for element in indices:
        atomsRemove.append(int(element))

    sorted(atomsRemove, reverse=True)

    del file[atomsRemove]

    return file


slab = fcc100('Pd', size=(4, 4, 5), periodic=True, vacuum=10)
slab = Remove(slab)
view(slab)