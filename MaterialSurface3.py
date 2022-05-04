########################################################################################################################
# Script to compute the Area of a slab. Uses a tilling scheme to compute area, taking in account non-planarity of a    #
# surface.                                                                                                             #
# Version 1.0 created on xx/xx/2021 by Alexandre Boucher. Last update: xx/xx/20xx. Version xx.xx                       #
########################################################################################################################
"""
So far, these are the materials supported by this program:
- cubic CeO2,
- MgO,
- CaNH (rocksalt)
- a-Quartz,
- any fcc metal (ensure they exist in AtomColor function),
- any bcc metal (Fe, W only so far, ad more by setting new atom type in AtomColor function),
- any hcp metal (same as fcc & bcc, ensure it exists in AtomColor).
Failures:
- hexagonal ZnO.

In order to compute surface energy, if needed, set getSurfaceEnergy = True line 389, and set appropriate files.

For materials' surface saturated in Hydrogen atoms, those atoms are removed in the tiling scheme, assuming their
very mobile nature won't provide a good enough approximation of the area measurement.
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from ase.visualize import view
from ase.io import read
from ase import Atoms, neighborlist
from ase.neighborlist import neighbor_list, natural_cutoffs
from ase.build import surface, bulk
from ase.spacegroup import get_spacegroup, crystal
from shapely.geometry import Polygon

matplotlib.use('TkAgg')
plt.ion()

"""
Set up all the path toward the required filed.
/home/alexbou/Tungstene_Esurf/W110_3x4x5_3vac/
/home/alexbou/Tungstene_Esurf/W110_3x4x5_0R/
/home/alexbou/Tungstene_Esurf/TungsteneDFTparameters/isif3/
"""

path_rel = '/home/alexbou/Ru_fcc_611/1x3/R3_L6/'  # path toward the relaxed slab.
path_froz = '/home/alexbou/Ru_fcc_611/1x3/R0_L6/'  # path toward the frozen slab.
path_bulk = '/home/alexbou/Bulk/fcc/'  # Path toward the bulk OUTCAR.

getSurfaceEnergy = True  # Compute surface energy or just surface?
testError = True  # Ask if error while generating initial tiles to remove them all.

# Read CONTCAR relaxed slab.
bulkmat = read(path_bulk + 'CONTCAR')
slab = read(path_rel + 'CONTCAR')
view(slab)

"""
Generates all the functions we need to run the script.
"""


def BuildDictionary(bulkmat):
    bulkcell = bulkmat.get_cell()
    bulkmat = Atoms([atom for atom in bulkmat], cell=bulkcell, pbc=True) * (2, 2, 2)

    # build a list of atoms present in the structure.
    bulkAtoms = []
    for atom in bulkmat:
        if atom.symbol not in bulkAtoms:
            bulkAtoms.append(atom.symbol)

    # Build all possible pair of atoms in the system.
    atomPairs = []
    for element1 in bulkAtoms:
        for element2 in bulkAtoms:
            pair = [element1, element2]
            if pair not in atomPairs:
                atomPairs.append(pair)

    # Build a cufOff for each atom pair in the system. Inject it into a dictionary.
    dictionary = {}
    for pair in atomPairs:
        d = 100
        for atom1 in bulkmat:
            if atom1.symbol == pair[0]:
                for atom2 in bulkmat:
                    if atom1.index != atom2.index and atom2.symbol == pair[1]:
                        distance = np.linalg.norm(atom2.position - atom1.position)
                        if distance < d:
                            d = distance
        dictionary.update({(pair[0], pair[1]): round(d * 1.09, 3)})

    return bulkAtoms, dictionary


def NaturalCutOff(bulkmat, atomList, dictionary):
    bulkcell = bulkmat.get_cell()
    bulkmat = Atoms([atom for atom in bulkmat], cell=bulkcell, pbc=True)

    neighList = neighbor_list('i', a=bulkmat, cutoff=dictionary, self_interaction=False)
    cnList = np.bincount(neighList)

    cnMax = []
    for element in atomList:
        for atom in bulkmat:
            if atom.symbol == element and cnList[atom.index] not in cnMax:
                cnMax.append(cnList[atom.index])

    return cnMax


def SurfaceAtoms(mat, elements, cutOffList, maxCN, side):
    """
    :param mat: The slab understudied.
    :param elements: list containing elements in the slab, provided by NaturalCutOff (above).
    :param maxCN: list containing bulk max CN for each element, provided by NaturalCutOff (above).
    :param side: top or bottom side of the slide.
    :return: List of atom either on top or on the bottom of the slab (undercoodinated).
    """
    # Reproduce cell
    cell = mat.get_cell()
    unit = Atoms([i for i in mat], cell=cell, pbc=True) * (2, 2, 1)

    # Get zMax.
    zMax = 0
    for atom in unit:
        if atom.position[2] > zMax:
            zMax = atom.position[2]

    surfaceAtoms = []
    i = neighborlist.neighbor_list('i', unit, cutOffList, self_interaction=False)
    CN = np.bincount(i)

    for atom in unit:
        if side == 'top':
            if atom.position[2] > (zMax / 2):
                atomIndex = atom.index
                atomCN = CN[atomIndex]
                atomSymbol = atom.symbol
                for i in elements:
                    if i == atomSymbol:
                        index = elements.index(atomSymbol)
                        if atomCN < maxCN[index]:
                            surfaceAtoms.append(atom)
        if side == 'bot':
            if atom.position[2] < (zMax / 2):
                atomIndex = atom.index
                atomCN = CN[atomIndex]
                atomSymbol = atom.symbol
                for i in elements:
                    if i == atomSymbol:
                        index = elements.index(atomSymbol)
                        if atomCN < maxCN[index]:
                            surfaceAtoms.append(atom)

    # Shift them so that smalled x-pos and y_pos is at coordinates (0, 0, z).
    xMin = 10e89
    yMin = 10e89
    for atom in surfaceAtoms:
        if atom.position[0] < xMin and atom.position[1] < yMin:
            xMin = atom.position[0]
            yMin = atom.position[0]
    for atom in surfaceAtoms:
        atom.position = atom.position - [xMin, yMin, 0]

    return surfaceAtoms


def DelHydrogen(surfAtoms, side):
    surfaceAtomsCorrected = []
    if side == 'top':
        for atom in surfAtoms:
            if atom.symbol != 'H':
                surfaceAtomsCorrected.append(atom)
    if side == 'bot':
        for atom in surfAtoms:
            if atom.symbol != 'H':
                surfaceAtomsCorrected.append(atom)

    return surfaceAtomsCorrected


def PlotCell(mat, surfAtoms, side):
    """
    :param mat: The slab we're working on.
    :param surfAtoms: The list of surface atoms, provided by SurfaceAtoms function (above).
    :param side: The side we working on (top or bot).
    :return: Plots the contour of the supercell.
    """
    # Get Zaverage.
    zAverage = sum([atom.position[2] for atom in surfAtoms]) / len(surfAtoms)

    # Get cell.
    cell = mat.get_cell()
    a = cell[0]  # First vector.
    b = cell[1]  # Second vector.
    c = cell[2]  # Third vector.

    # For top side it can be a bit more complex if cell[2] has x and z component.
    if side == 'top':
        zMax = 0
        for atom in mat:
            if atom.position[2] > zMax:
                zMax = atom.position[2]
        # Take the ratio zMax/(Znorm of cell[2]). Multiply ratio by c[0] and c[1] to get shifting.
        shiftRatio = zMax / c[2]
        xShift = shiftRatio * c[0]
        yShift = shiftRatio * c[1]

    # For bot side its easier.
    if side == 'bot':
        xShift = 0
        yShift = 0

    ax.quiver(0 + xShift, 0 + yShift, zAverage,
              a[0] + xShift, a[1] + yShift, 0, length=np.linalg.norm(a), arrow_length_ratio=0.00, color="b",
              lw=3, alpha=0.5, normalize=True)
    ax.quiver(0 + xShift, 0 + yShift, zAverage,
              b[0] + xShift, b[1] + yShift, 0, length=np.linalg.norm(b), arrow_length_ratio=0.00, color="b",
              lw=3, alpha=0.5, normalize=True)
    ax.quiver(a[0] + xShift, a[1] + yShift, zAverage,
              b[0] + xShift, b[1] + yShift, 0, length=np.linalg.norm(b), arrow_length_ratio=0.00, color="b",
              lw=3, alpha=0.5, normalize=True)
    ax.quiver(b[0] + xShift, b[1] + yShift, zAverage,
              a[0] + xShift, a[1] + yShift, 0, length=np.linalg.norm(a), arrow_length_ratio=0.00, color="b",
              lw=3, alpha=0.5, normalize=True)

    return


# Provide size and color for each atom type needed.
def AtomColor(atom):
    if atom.symbol == 'Si':
        return '#FF9900', 150, 0.4
    if atom.symbol == 'O':
        return '#FF0000', 100, 0.4
    if atom.symbol == 'Mg':
        return '#FF6600', 200, 0.4
    if atom.symbol == 'Ce':
        return '#FFFF99', 210, 0.4
    if atom.symbol == 'Zn':
        return '#6666FF', 170, 0.4
    if atom.symbol == 'Pd':
        return '#9999CC', 250, 0.4
    if atom.symbol == 'Cu':
        return '#FF6633', 250, 0.4
    if atom.symbol == 'Fe':
        return '#666666', 250, 0.4
    if atom.symbol == 'Ni':
        return '#999999', 250, 0.4
    if atom.symbol == 'Au':
        return '#ffcc00', 250, 0.4
    if atom.symbol == 'Ru':
        return '#3399CC', 250, 0.4
    if atom.symbol == 'Ni':
        return '#999999', 250, 0.4
    if atom.symbol == 'Co':
        return '#FF99CC', 250, 0.4
    if atom.symbol == 'W':
        return '#33FFFF', 250, 0.4
    if atom.symbol == 'Ca':
        return '#00ff00', 450, 0.4
    if atom.symbol == 'N':
        return '#0033ff', 90, 0.4
    if atom.symbol == 'H':
        return '#cccccc', 50, 0.4


def TriangleArea(x, y, z):
    V1 = x - y
    D1 = np.linalg.norm(V1)
    V2 = y - z
    D2 = np.linalg.norm(V2)
    V3 = x - z
    D3 = np.linalg.norm(V3)
    SP = (D1 + D2 + D3) / 2
    Area = (SP * (SP - D1) * (SP - D2) * (SP - D3)) ** 0.5
    return Area


# Eliminates extra atoms from the (2, 2, 1) mult (line 53).
def Control(mat, surfAtoms, side):
    # Get zAverage.
    zAverage = sum([atom.position[2] for atom in surfAtoms]) / len(surfAtoms)

    # Check flat area. ControlArea will be used to discriminate surface atoms.
    cell = mat.get_cell()
    flatMatrix = np.array([[[cell[0][0], cell[0][1]], [cell[1][0], cell[1][1]]]])
    controlArea = np.linalg.det(flatMatrix)

    a = cell[0]  # First vector.
    b = cell[1]  # Second vector.
    c = cell[2]  # Third vector.

    # Once again, more tricky of side = top and cell[2] has x and y components!
    if side == 'top':
        zMax = 0
        for atom in mat:
            if atom.position[2] > zMax:
                zMax = atom.position[2]
        # Take the ratio zMax/(Znorm of cell[2]). Multiply ratio by c[0] and c[1] to get shifting.
        shiftRatio = zMax / c[2]
        xShift = shiftRatio * c[0]
        yShift = shiftRatio * c[1]

    # For bot side its easier.
    if side == 'bot':
        xShift = 0
        yShift = 0

    Cell_corner1 = np.array([0 + xShift, 0 + yShift, zAverage])
    Cell_corner2 = np.array([cell[0][0] + xShift, cell[0][1] + yShift, zAverage])
    Cell_corner3 = np.array([cell[1][0] + xShift, cell[1][1] + yShift, zAverage])
    Cell_corner4 = np.array([cell[0][0] + cell[1][0] + xShift, cell[0][1] + cell[1][1] + yShift, zAverage])

    cleanSurfAtoms = []
    for atom in surfAtoms:
        position = atom.position
        Area1 = TriangleArea(position, Cell_corner1, Cell_corner2)
        Area2 = TriangleArea(position, Cell_corner1, Cell_corner3)
        Area3 = TriangleArea(position, Cell_corner2, Cell_corner4)
        Area4 = TriangleArea(position, Cell_corner3, Cell_corner4)
        AreaTest = Area1 + Area2 + Area3 + Area4
        if AreaTest <= 1.3 * controlArea:
            cleanSurfAtoms.append(atom)

    return cleanSurfAtoms


def PlotAtoms(atomList):
    for atom in atomList:
        color, size, alpha = AtomColor(atom)
        ax.scatter3D(atom.position[0], atom.position[1], atom.position[2], c=color, s=size, alpha=alpha, edgecolors='k')
    return


def PseudoBonds(atomList, elements):
    cutOffs = natural_cutoffs(atoms=atomList, mult=1.7)
    # Generate pseudo bonds
    for atom1 in atomList:
        atomSymbol = atom1.symbol
        index = elements.index(atomSymbol)
        cutOff = cutOffs[index] * 1.5  # Allow a slightly bigger cut off - its only pseudo bonds!
        for atom2 in atomList:
            if atom1.index != atom2.index:
                vector = atom2.position - atom1.position
                distance = np.linalg.norm(vector)
                if distance <= cutOff:
                    ax.quiver(atom1.position[0], atom1.position[1], atom1.position[2],
                              vector[0], vector[1], vector[2], length=1, arrow_length_ratio=0.09,
                              color='k', lw=2, normalize=True, alpha=0.3)
    return


def PlotLabel(atomList):
    for atom in atomList:
        label = str(atom.index)
        ax.text(atom.position[0] + (0.0009 * atom.position[0]),
                atom.position[1] + (0.0009 * atom.position[1]),
                atom.position[2] + (0.0009 * atom.position[2]),
                label,
                c="k",
                alpha=1)
    return


def PreAddedTiles(elements, surfaceAtoms):
    cutOff = natural_cutoffs(atoms=surfaceAtoms, mult=1.3)
    cutOff = (sum(i for i in cutOff) / len(cutOff)) * 1.9
    tilesTotal = []
    for i in range(0, len(elements)):
        element = elements[i]

        for atomI in surfaceAtoms:
            if atomI.symbol == element:
                positionI = atomI.position
                for atomJ in surfaceAtoms:
                    if atomI.index != atomJ.index:
                        positionJ = atomJ.position
                        vectorIJ = positionJ - positionI
                        distanceIJ = np.linalg.norm(vectorIJ)
                        if distanceIJ <= cutOff:
                            for atomK in surfaceAtoms:
                                if atomK.index != atomI.index and atomK.index != atomJ.index:
                                    positionK = atomK.position
                                    vectorIK = positionK - positionI
                                    vectorJK = positionK - positionJ
                                    distanceIK = np.linalg.norm(vectorIK)
                                    distanceJK = np.linalg.norm(vectorJK)
                                    if distanceJK <= cutOff and distanceIK <= cutOff:
                                        tile = [atomI.index, atomJ.index, atomK.index]
                                        tile.sort()
                                        if tile not in tilesTotal:
                                            tilesTotal.append(tile)
    # Test.
    for tile1 in tilesTotal:
        positionref = []
        for i in tile1:
            for atom in surfaceAtoms:
                if atom.index == i:
                    positionref.append(atom.position)

        for tile2 in tilesTotal:
            if not np.array_equal(tile1, tile2):
                positionTest = []
                for j in tile2:
                    for atom in surfaceAtoms:
                        if atom.index == j:
                            positionTest.append(atom.position)

                p1 = Polygon([positionref[0], positionref[1], positionref[2]])
                p2 = Polygon([positionTest[0], positionTest[1], positionTest[2]])
                if p1.overlaps(p2):
                    tilesTotal.remove(tile2)
    return tilesTotal


def AddTile(tile, surfaceAtoms):
    atomsPositions = []
    for i in tile:
        for atom in surfaceAtoms:
            if atom.index == i:
                atomsPositions.append(atom.position)
    # Coordinates.
    X = np.array([atomsPositions[0][0], atomsPositions[1][0], atomsPositions[2][0]])
    Y = np.array([atomsPositions[0][1], atomsPositions[1][1], atomsPositions[2][1]])
    Z = np.array([atomsPositions[0][2], atomsPositions[1][2], atomsPositions[2][2]])

    # Z_average is ised to define the color of the tile.
    Zweight = [atomsPositions[0][2], atomsPositions[1][2], atomsPositions[2][2]]
    Zweight.sort()
    zmax = max(Zweight)
    zmin = min(Zweight)
    if (zmax - zmin) > 0.5:
        Zaverage = (sum(Zweight) / 3) / (2 * (zmax - zmin))
    else:
        Zaverage = (sum(Zweight) / 3)

    verts = [list(zip(X, Y, Z))]

    ax.add_collection3d(
        Poly3DCollection(verts, edgecolors="k", lw=0.1, facecolors=plt.cm.viridis(Zaverage), alpha=0.2),
        zdir="z")
    return


########################################################################################################################
# Compute surface.                                                                                                     #
########################################################################################################################


# The different supported materials.
# MgO: https://materialsproject.org/materials/mp-1265/
MgO = read('./CIF_&_CONTCAR_files/MgO.cif')
# alpha-Quartz: https://materialsproject.org/materials/mp-6930/
aQuartz = read('./CIF_&_CONTCAR_files/SiO2_mp-6930_primitive.cif')
# cubic ceria: https://materialsproject.org/materials/mp-20194/
cubicCeria = read('./CIF_&_CONTCAR_files/mp-20194_CeO2_cubic.cif')
# hexagonal ZnO: https://materialsproject.org/materials/mp-2133/
hexZnO = read('./CIF_&_CONTCAR_files/mp-2133_ZnO_hex.cif')
# Pd, or any fcc metal.
#fccMetal = bulk('Cu')
# bcc metals, Iron.
a = 2.87
Fe = crystal('Fe', [(0, 0, 0)], spacegroup=229, cellpar=[a, a, a, 90, 90, 90])
# Read a surface on local laptop.

########################################################################################################################
# Set files needed.                                                                                                    #
########################################################################################################################

if getSurfaceEnergy:
    matOUTCARrelax = read(path_rel + 'OUTCAR')  # OUTCAR relaxed slab.
    matEnergyRelax = matOUTCARrelax.get_potential_energy()
    matOUTCARfrozen = read(path_froz + 'OUTCAR')  # Path toward frozen OUTCAR file.
    matEnergyFrozen = matOUTCARfrozen.get_potential_energy()
    matBulk = read(path_bulk + 'OUTCAR')  # Path toward bulk OUTCAR.
    nBulk = 0
    for atom in matBulk:
        nBulk += 1
    matEnergyBulk = matBulk.get_potential_energy() / nBulk

########################################################################################################################
#                                                                                                                      #
########################################################################################################################
# Step1: Get a list of atom's elements and a dictionary of atomic distances.
elements, dictionary = BuildDictionary(bulkmat=bulkmat)

# Step2: Get CN max of each atom in the system.
maxCN = NaturalCutOff(bulkmat=bulkmat, atomList=elements, dictionary=dictionary)

# Get top and bottom atoms.
topAtoms = SurfaceAtoms(mat=slab, elements=elements, cutOffList=dictionary, maxCN=maxCN, side='top')
topAtoms = DelHydrogen(surfAtoms=topAtoms, side='top')
botAtoms = SurfaceAtoms(mat=slab, elements=elements, cutOffList=dictionary, maxCN=maxCN, side='bot')
botAtoms = DelHydrogen(surfAtoms=botAtoms, side='bot')

# Remove extra atoms.
cleanTopAtoms = Control(mat=slab, surfAtoms=topAtoms, side='top')
cleanBotAtoms = Control(mat=slab, surfAtoms=botAtoms, side='bot')

# Generate intitial tiles.
initialTilesTop = PreAddedTiles(elements, cleanTopAtoms)
initialTilesBot = PreAddedTiles(elements, cleanBotAtoms)

########################################################################################################################
# Top surface business.                                                                                                #
########################################################################################################################
fig = plt.figure(figsize=(10, 10), clear=True)
ax = fig.add_subplot(111, projection='3d')
plt.axis(False)

# Top cell contour and top atoms.
PlotAtoms(cleanTopAtoms)
PseudoBonds(atomList=cleanTopAtoms, elements=elements)
PlotCell(mat=slab, surfAtoms=topAtoms, side='top')
PlotLabel(cleanTopAtoms)

# Plot initial tiles.
for tile in initialTilesTop:
    AddTile(tile, cleanTopAtoms)

zAverage = sum([atom.position[2] for atom in cleanTopAtoms]) / len(cleanTopAtoms)
ax.set_zlim3d(zAverage - 5, zAverage + 5)
plt.show()

'''
Now, add Tiles if needed.
'''
errorTiling = input('Is there an error in tile generation? (y/n)\n')
if errorTiling == 'y':
    initialTilesTop = []
    ax.clear()
    plt.axis(False)

    # Top cell contour and top atoms.
    PlotAtoms(cleanTopAtoms)
    PseudoBonds(atomList=cleanTopAtoms, elements=elements)
    PlotCell(mat=slab, surfAtoms=topAtoms, side='top')
    PlotLabel(cleanTopAtoms)

    # Plot initial tiles.
    for tile in initialTilesTop:
        AddTile(tile, cleanTopAtoms)

    zAverage = sum([atom.position[2] for atom in cleanTopAtoms]) / len(cleanTopAtoms)
    ax.set_zlim3d(zAverage - 5, zAverage + 5)

    plt.show()

totTilesTop = initialTilesTop
addTopTile = 'y'
while addTopTile == 'y':
    addTopTile = input('Would you like to cover any other area (y/n)?\n')
    if addTopTile == 'y':
        TILE = input('Enter indices of corners (a, b, c):\n').split()
        tile = []
        for element in TILE:
            for atom in cleanTopAtoms:
                if str(atom.index) == element:
                    tile.append(atom.index)
        tile.sort()
        AddTile(tile, cleanTopAtoms)
        totTilesTop.append(tile)
    if len(addTopTile) > 2:
        TILE = addTopTile.split()
        tile = []
        for element in TILE:
            for atom in cleanTopAtoms:
                if str(atom.index) == element:
                    tile.append(atom.index)
        tile.sort()
        AddTile(tile, cleanTopAtoms)
        totTilesTop.append(tile)
        addTopTile = 'y'

'''
Do you need to remove Tiles?
'''
for tile in totTilesTop:
    tile.sort()
removeTopTile = 'y'
while removeTopTile == 'y':
    removeTopTile = input('Would you like to remove any area (y/n)?\n')
    if removeTopTile == 'y':
        removeTopTile = input('Enter indices of corners (a, b, c):\n')
    if len(removeTopTile) > 2:
        removeTopTile = removeTopTile.split()
        removeTile = []

        for element in removeTopTile:
            for atom in cleanTopAtoms:
                if str(atom.index) == element:
                    removeTile.append(atom.index)
        removeTile.sort()
        for tile in totTilesTop:
            if removeTile == tile:
                totTilesTop.remove(tile)

        ax.clear()
        plt.axis(False)

        # Top cell contour and top atoms.
        PlotAtoms(cleanTopAtoms)
        PseudoBonds(atomList=cleanTopAtoms, elements=elements)
        PlotCell(mat=slab, surfAtoms=topAtoms, side='top')
        PlotLabel(cleanTopAtoms)

        # Plot initial tiles.
        for tile in totTilesTop:
            AddTile(tile, cleanTopAtoms)

        zAverage = sum([atom.position[2] for atom in cleanTopAtoms]) / len(cleanTopAtoms)
        ax.set_zlim3d(zAverage - 5, zAverage + 5)

        plt.show()
        removeTopTile = 'y'

'''
Just in case.
'''
totTilesTop = initialTilesTop
addTopTile = 'y'
while addTopTile == 'y':
    addTopTile = input('Would you like to cover any other area (y/n)?\n')
    if addTopTile == 'y':
        TILE = input('Enter indices of corners (a, b, c):\n').split()
        tile = []
        for element in TILE:
            for atom in cleanTopAtoms:
                if str(atom.index) == element:
                    tile.append(atom.index)
        tile.sort()
        AddTile(tile, cleanTopAtoms)
        totTilesTop.append(tile)
    if len(addTopTile) > 2:
        TILE = addTopTile.split()
        tile = []
        for element in TILE:
            for atom in cleanTopAtoms:
                if str(atom.index) == element:
                    tile.append(atom.index)
        tile.sort()
        AddTile(tile, cleanTopAtoms)
        totTilesTop.append(tile)
        addTopTile = 'y'

########################################################################################################################
# Bot surface business.                                                                                                #
########################################################################################################################
fig = plt.figure(figsize=(10, 10), clear=True)
ax = fig.add_subplot(1, 1, 1, projection='3d')
plt.axis(False)

# Bot cell contour and bot atoms.
PlotCell(mat=slab, surfAtoms=botAtoms, side='bot')
PlotAtoms(cleanBotAtoms)
PseudoBonds(atomList=cleanBotAtoms, elements=elements)
PlotLabel(cleanBotAtoms)

# Plot initial tiles.
for tile in initialTilesBot:
    AddTile(tile, cleanBotAtoms)

zAverage = sum([atom.position[2] for atom in cleanBotAtoms]) / len(cleanBotAtoms)
ax.set_zlim3d(zAverage - 5, zAverage + 5)
plt.show()

'''
Now, add Tiles if needed.
'''
errorTiling = input('Is there an error in tile generation? (y/n)\n')
if errorTiling == 'y':
    initialTilesBot = []
    ax.clear()
    plt.axis(False)

    # Top cell contour and top atoms.
    PlotAtoms(cleanBotAtoms)
    PseudoBonds(atomList=cleanBotAtoms, elements=elements)
    PlotCell(mat=slab, surfAtoms=botAtoms, side='top')
    PlotLabel(cleanBotAtoms)

    # Plot initial tiles.
    for tile in initialTilesBot:
        AddTile(tile, cleanBotAtoms)

    zAverage = sum([atom.position[2] for atom in cleanBotAtoms]) / len(cleanBotAtoms)
    ax.set_zlim3d(zAverage - 5, zAverage + 5)

    plt.show()

totTilesBot = initialTilesBot
addBotTile = 'y'
while addBotTile == 'y':
    addBotTile = input('Would you like to cover any other area (y/n)?\n')
    if addBotTile == 'y':
        TILE = input('Enter indices of corners (a, b, c):\n').split()
        tile = []
        for element in TILE:
            for atom in cleanBotAtoms:
                if str(atom.index) == element:
                    tile.append(atom.index)
        tile.sort()
        print(tile)
        AddTile(tile, cleanBotAtoms)
        totTilesBot.append(tile)
    if len(addBotTile) > 2:
        TILE = addBotTile.split()
        tile = []
        for element in TILE:
            for atom in cleanBotAtoms:
                if str(atom.index) == element:
                    tile.append(atom.index)
        tile.sort()
        AddTile(tile, cleanBotAtoms)
        totTilesBot.append(tile)
        addBotTile = 'y'

'''
Do you need to remove Tiles?
'''
for tile in totTilesBot:
    tile.sort()
removeBotTile = 'y'
while removeBotTile == 'y':
    removeBotTile = input('Would you like to remove any area (y/n)?\n')
    if removeBotTile == 'y':
        removeBotTile = input('Enter indices of corners (a, b, c):\n')
    if len(removeBotTile) > 2:
        removeBotTile = removeBotTile.split()
        removeTile = []

        for element in removeBotTile:
            for atom in cleanBotAtoms:
                if str(atom.index) == element:
                    removeTile.append(atom.index)
        removeTile.sort()
        for tile in totTilesBot:
            if removeTile == tile:
                totTilesBot.remove(tile)

        ax.clear()
        plt.axis(False)

        # Top cell contour and top atoms.
        PlotAtoms(cleanBotAtoms)
        PseudoBonds(atomList=cleanBotAtoms, elements=elements)
        PlotCell(mat=slab, surfAtoms=botAtoms, side='top')
        PlotLabel(cleanBotAtoms)

        # Plot initial tiles.
        for tile in totTilesBot:
            AddTile(tile, cleanBotAtoms)

        zAverage = sum([atom.position[2] for atom in cleanBotAtoms]) / len(cleanBotAtoms)
        ax.set_zlim3d(zAverage - 5, zAverage + 5)

        plt.show()
        removeBotTile = 'y'

'''
Just in case.
'''
totTilesBot = initialTilesBot
addBotTile = 'y'
while addBotTile == 'y':
    addBotTile = input('Would you like to cover any other area (y/n)?\n')
    if addBotTile == 'y':
        TILE = input('Enter indices of corners (a, b, c):\n').split()
        tile = []
        for element in TILE:
            for atom in cleanBotAtoms:
                if str(atom.index) == element:
                    tile.append(atom.index)
        tile.sort()
        print(tile)
        AddTile(tile, cleanBotAtoms)
        totTilesBot.append(tile)
    if len(addBotTile) > 2:
        TILE = addBotTile.split()
        tile = []
        for element in TILE:
            for atom in cleanBotAtoms:
                if str(atom.index) == element:
                    tile.append(atom.index)
        tile.sort()
        AddTile(tile, cleanBotAtoms)
        totTilesBot.append(tile)
        addBotTile = 'y'

########################################################################################################################
# Now, calculate surface.                                                                                              #
########################################################################################################################
# First top surface.
surfTop = 0
for tile in totTilesTop:
    positions = []
    for atom in cleanTopAtoms:
        if atom.index == tile[0]:
            positions.append(atom.position)
        if atom.index == tile[1]:
            positions.append(atom.position)
        if atom.index == tile[2]:
            positions.append(atom.position)

    surfTop += TriangleArea(positions[0], positions[1], positions[2])
print('Total top surface: ', surfTop, 'A²')

# Then, bot surface.
surfBot = 0
for tile in totTilesBot:
    positions = []
    for atom in cleanBotAtoms:
        if atom.index == tile[0]:
            positions.append(atom.position)
        if atom.index == tile[1]:
            positions.append(atom.position)
        if atom.index == tile[2]:
            positions.append(atom.position)

    surfBot += TriangleArea(positions[0], positions[1], positions[2])
print('Total bot surface: ', surfBot, 'A²')

########################################################################################################################
# Get surface energy.                                                                                                  #
########################################################################################################################
if getSurfaceEnergy:
    # Number of atoms in relaxed structure.
    nRelaxed = 0
    for atom in matOUTCARrelax:
        nRelaxed += 1
    # Number of atoms in frozen structure.
    nFrozen = 0
    for atom in matOUTCARfrozen:
        nFrozen += 1
    # Get Frozen surface energy.
    yFrozen = (matEnergyFrozen - (nFrozen * matEnergyBulk)) / (2 * surfBot)
    yRelaxed = (matEnergyRelax - (nRelaxed * matEnergyBulk) - (surfBot * yFrozen)) / surfTop
    # Convert from eV/A² to J/m²
    yFrozen = yFrozen * 1.602180000067e-19 * 10e+19
    yRelaxed = yRelaxed * 1.602180000067e-19 * 10e+19
    print('Bottom surface energy, J/m² = ', yFrozen)
    print('Top surface energy, J/m² =', yRelaxed)