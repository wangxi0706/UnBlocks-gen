# UnBlocks-Gen: 3D rock mass generator and analyser
# Copyright (C) 2020  Leandro Lima Rasmussen
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from unblocks import *

# dfn1 = DFN()
# dfn1.set_FirstRecBlock([0, 0, 0], [5, 5, 1.])
# dfn1.add_FractureSet()

# dfn2=DFN()
# dfn2.set_FirstRecBlock([0,0,1],[5,5,3])
# dfn2.add_FractureSet()
# dfn2.fractureSets[0].add_CircularFracture([2.5, 2.5, 2], 0, 0, 1)


# generator = Generator()
# generator.add_Fixed_Region([-0.01, -0.01, -0.01], [5.1, 5.1, 1.01])
# generator.generate_RockMass_Multi(dfn1)
# generator.generate_RockMass_Multi(dfn2)

dfn=DFN()
n=11
dfn.set_RegionMaxCorner([8,8,n])
dfn.add_FractureSet()
generator = Generator()
generator.add_Fixed_Region([-0.01, -0.01, -0.01], [10.1, 10.1, 1.01])
for i in range(n-1):
    dfn.fractureSets[0].add_CircularFracture([5, 5, i+1], 0, 0, 1)
generator.generate_RockMass_Multi(dfn)

s="DataSaveGeo_BlkCut_"+str(n)
generator.export_BlocksDDAOpt(s)
generator.export_BlocksVtk(s)

# for i in range(len(generator.blocks)):
# generator.blocks[i].export_BlockVtk("Block"+str(i))
