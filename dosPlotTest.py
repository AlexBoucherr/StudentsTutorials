import pyprocar
import matplotlib.pyplot as plt
from ase.io import read
from ase.io.vasp import read_vasp_xml
import os

run = '/home/alex/pdh_oho_bridge_singlepoint/vasprun.xml'

# 22, 29, 28, 43, 18, 56, 68 Oxygens
# 97 Palladium.
# 73, 71
fig, ax = pyprocar.dosplot(filename=run,
                           mode='stack_orbitals',
                           atoms=[28, 97],
                           elimit=[-7, 2],
                           dos_limit=[-10, 10],
                           orientation='horizontal',
                           plot_total=False,
                           title=r'Title')
fig.show()
