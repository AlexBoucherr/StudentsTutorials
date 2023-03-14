from ase.io import read
import shutil

POSCAR = read('./POSCAR')

Elements = []
for i in POSCAR:
    if i.symbol not in Elements:
        Elements.append(i.symbol)

POTCARS_needed = []

NAME = 'POTCAR_'

for i in Elements:
    j = 'POTCAR_' + i
    NAME = NAME + str(i)
    POTCARS_needed.append(j)

Path = './' + NAME

with open(Path, 'w') as outfile:
    for i in POTCARS_needed:
        path2 = '/home/c.c1981790/POTCAR_old/' + i
        with open(path2) as infile:
            for line in infile:
                outfile.write(line)

shutil.move(NAME, './POTCAR_gen/')
