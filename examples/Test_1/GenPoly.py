from unblocks import *
import generateBlock as GB
import numpy as np

################################
# For block B
w = 30  # x
l = 10  # y
h = 20  # z

tol = 1
w += tol
l += tol
h += tol

vertexes = []
vertexes.append([tol, l+tol, h+tol])
vertexes.append([tol, tol, h+tol])
vertexes.append([w+tol, tol, h+tol])
vertexes.append([w+tol, l+tol, h+tol])
vertexes.append([tol, l+tol, tol])
vertexes.append([tol, tol, tol])
vertexes.append([w+tol, tol, tol])
vertexes.append([w+tol, l+tol, tol])
edges = [[0, 1], [1, 2], [2, 3], [3, 0], [0, 4], [1, 5], [2, 6], [3, 7],
         [4, 5], [5, 6], [6, 7], [7, 4]]
faces = [[0, 1, 2, 3], [0, 4, 5, 1], [0, 3, 7, 4], [2, 6, 7, 3],
         [1, 5, 6, 2], [4, 7, 6, 5]]
block1 = GB.Block(vertexes, edges, faces, fixed=1)

GB.genBlkTet([block1], "cuboid.poly")

# get the fixed boundary, 1 bottom, 2 left, 3 front, 4 behind
# left block
generator = Generator()
dfn = DFN()
min = [0, 0, 0]
max = [0+tol, l+2*tol, h+tol]
dfn.set_FirstRecBlock(min, max)
dfn.add_FractureSet()
generator.generate_RockMass_Multi(dfn)
# bottom block
dfn = DFN()
min = [tol, 0, 0]
max = [4*w+tol, l+2*tol, tol]
dfn.set_FirstRecBlock(min, max)
dfn.add_FractureSet()
generator.generate_RockMass_Multi(dfn)
# front block
dfn = DFN()
min = [tol, 0, tol]
max = [4*w+tol, tol, h+tol]
dfn.set_FirstRecBlock(min, max)
dfn.add_FractureSet()
generator.generate_RockMass_Multi(dfn)
# back block
dfn = DFN()
min = [tol, l+tol, tol]
max = [4*w+tol, l+2*tol, h+tol]
dfn.set_FirstRecBlock(min, max)
dfn.add_FractureSet()
generator.generate_RockMass_Multi(dfn)

generator.add_Fixed_Region([-100*w, -100*w, -100*w], [100*w, 100*w, 100*w])
#################################################
s = "FixBnd"
generator.export_BlocksVtk(s)
generator.export_BlocksDDAOpt(s)
