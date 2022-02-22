from unblocks import *
import generateBlock as GB
import numpy as np

################################
# For block B
w = 3  # x
l = 3  # y
h = 3  # z

tol = 0.5
w += tol
l += tol
h += tol

vertexes = []
vertexes.append([0.5*w+tol, 0.5*l+l+tol, 0.1+h+tol])
vertexes.append([0.5*w+tol, 0.5*l+tol, 0.1+h+tol])
vertexes.append([0.5*w+w+tol, 0.5*l+tol, 0.1+h+tol])
vertexes.append([0.5*w+w+tol, 0.5*l+l+tol, 0.1+h+tol])
vertexes.append([0.5*w+tol, 0.5*l+l+tol, 0.1+tol])
vertexes.append([0.5*w+tol, 0.5*l+tol, 0.1+tol])
vertexes.append([0.5*w+w+tol, 0.5*l+tol, 0.1+tol])
vertexes.append([0.5*w+w+tol, 0.5*l+l+tol, 0.1+tol])
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
# generator.generate_RockMass_Multi(dfn)
# bottom block
dfn = DFN()
min = [tol, 0, 0]
max = [2*w+tol, 2*l+2*tol, tol]
dfn.set_FirstRecBlock(min, max)
dfn.add_FractureSet()
generator.generate_RockMass_Multi(dfn)
# front block
dfn = DFN()
min = [tol, 0, tol]
max = [4*w+tol+0.2, tol, h+tol+0.2]
dfn.set_FirstRecBlock(min, max)
dfn.add_FractureSet()
# generator.generate_RockMass_Multi(dfn)
# back block
dfn = DFN()
min = [tol, l+tol, tol]
max = [4*w+tol, l+2*tol, h+tol]
dfn.set_FirstRecBlock(min, max)
dfn.add_FractureSet()
# generator.generate_RockMass_Multi(dfn)

generator.add_Fixed_Region([-100*w, -100*w, -100*w], [100*w, 100*w, 100*w])
#################################################
s = "FixBnd"
generator.export_BlocksVtk(s)
generator.export_BlocksDDAOpt(s)
