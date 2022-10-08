from unblocks import *
import generateBlock as GB
import numpy as np

################################
# For block B
w = 3  # x
l = 3  # y
h = 3  # z

generator = Generator()
dfn = DFN()

dfn.set_RegionMaxMinCorner(max, min)

# for i in range(7):
dfn.add_FractureSet()
#     print(i)
#     dfn.fractureSets[0].add_CircularFracture(
#         # [0.5*w+tol*1.1+np.tan(np.pi/4)*0.5*i, 0.5*l+tol*1.1, 0.1+tol*1.1], 0, 45, 30)
#         # [0.1+tol*1.1+np.tan(np.pi/4)*1.1*i, 0.1+tol*1.1, 0.1+tol*1.1], 0, 45, 1)
#         [min[0]+0.2+i*0.5, min[1]+0.1, min[2] + 0.2+i*0.5], 0, 45, 1)
#     # dfn.fractureSets[0].add_CircularFracture(
#     #     # [0.5*w+tol*1.1+np.tan(np.pi/4)*0.5*i, 0.5*l+tol*1.1, 0.1+tol*1.1], 0, 45, 30)
#     #     # [0.1+tol*1.1+np.tan(np.pi/4)*1.1*i, 0.1+tol*1.1, 0.1+tol*1.1], 0, 45, 1)
#     #     [w-0.2-i*0.5, 0.1, 0.1+i*0.5], 60, 45, 10)

generator.generate_RockMass(dfn)

# # bottom block
# dfn = DFN()
# min = [tol, 0, 0]
# max = [2*w+tol, 2*l+2*tol, tol]
# dfn.set_FirstRecBlock(min, max)
# dfn.add_FractureSet()
# generator.generate_RockMass_Multi(dfn)

# generator.add_Fixed_Region(
#     [tol-0.1, 0-0.1, 0-0.1], [2*w+tol+0.1, 2*l+2*tol+0.1, tol+0.1])
#################################################
s = "Polyhedra"
generator.export_BlocksVtk(s)
generator.export_BlocksDDAOpt(s)
