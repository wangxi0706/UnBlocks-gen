
from unblocks import *
generator = Generator()

dfn = DFN()
dfn.set_FirstRecBlock([0,0,0], [6,6,1])
dfn.add_FractureSet()
generator.generate_RockMass_Multi(dfn)
generator.input_InputBound(0, [-1, 0, 0])  # 0=disp, 1=load

dfn = DFN()
dfn.set_FirstRecBlock([2, 2, 1], [4, 4, 2])
dfn.add_FractureSet()
generator.generate_RockMass_Multi(dfn)


generator.add_Fixed_Region([-.1, -.1, -.1], [6.1, 6.1, 1.1])
generator.export_BlocksDDAOpt("TestDynLoad")
generator.export_BlocksVtk("TestDynLoad")
