
from unblocks import *

dfn = DFN()
dfn.set_RegionMaxCorner([100, 100, 100])
dfn.add_FractureSet()
dfn.fractureSets[0].add_TriangularFracture(
    [0, 0, 50], [100, 0, 50], [0, 100, 50])

generator = Generator()
generator.generate_RockMass(dfn)

generator.export_BlocksDDAOpt("Simple_Cut")
generator.export_BlocksVtk("Simple_Cut")
