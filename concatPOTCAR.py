from ase.io import read

poscar = './MaterialGenerator/POSCAR_top_follow'  # The path toward the POSCAR file needing POTCAR.
potpath = './MaterialGenerator/POTCAR_old/'  # The path toward folder where all POTCARs are stored.

poscar = read(poscar)

# Check atoms present in poscar file.
atomType = []
for atom in poscar:
    if atom.symbol not in atomType:
        atomType.append(atom.symbol)

# Create a file called POTCAR.
potcar = open('POTCAR', 'w')

# Fill in POTCAR with all atomic POTCARs required for the calculation.
for element in atomType:
    potsearch = 'POTCAR_' + element
    potdata = potpath + potsearch
    with open(potdata, 'r') as singlepot:
        lines = singlepot.readlines()
        for line in lines:
            potcar.write(line)
        potcar.write('\n')

potcar.close()
