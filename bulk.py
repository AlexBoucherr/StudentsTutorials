from ase.build import bulk
from ase.visualize import view

"""
Use the bulk object to generate bulk material.
"""
# Bulk Palladium.
material = bulk('Pd', crystalstructure='fcc', cubic=False)
view(material)
