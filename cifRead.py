from ase.io import read
from ase import Atoms
from ase.visualize import view

"""
Read downloaded cif file, or simply create a cell yourself, using 'Atoms'.
"""

mgo = read('/home/alexbou/MgO.cif')
view(mgo)

graphene = Atoms('2C', scaled_positions=[[0.33, 0.66, 0.5], [0.66, 0.33, 0.5]], cell=(2.468, 2.468, 10.000, 90.0, 90.0, 120.0))
view(graphene)

graphene = graphene * (3, 3, 1)
view(graphene)