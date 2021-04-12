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
# dfn1.set_FirstRecBlock([0, 0, 0.1], [10, 10, 2.1])  # must be positive
# dfn1.add_FractureSet()
# dfn2 = DFN()
# dfn2.set_FirstRecBlock([2, 2, 2], [8, 8, 4])
# dfn2.add_FractureSet()

dfn1 = DFN()
dfn1.set_FirstRecBlock([-5, -5, 0], [5, 5, 2.])
dfn1.add_FractureSet()

dfn2 = DFN()
dfn2.set_FirstRecBlock([-5, 0, 2.], [5, 5, 4])
dfn2.add_FractureSet()
# dfn2.fractureSets[0].add_CircularFracture([0,0,3],0,90,0.5)

dfn3 = DFN()
dfn3.set_FirstRecBlock([-5, -5, 2.], [5, 0, 4])
dfn3.add_FractureSet()

dfn4 = DFN()
dfn4.set_FirstRecBlock([-5, -5, 4.], [5, 0, 6])
dfn4.add_FractureSet()

dfn5 = DFN()
dfn5.set_FirstRecBlock([-5, 0, 4.], [5, 5, 6])
dfn5.add_FractureSet()

generator = Generator()
generator.add_Fixed_Region([-5.1, -5.1, -0.1], [5.1, 5.1, 2.02])
generator.generate_RockMass_Multi(dfn1)
generator.generate_RockMass_Multi(dfn2)
generator.generate_RockMass_Multi(dfn3)
generator.generate_RockMass_Multi(dfn4)
# generator.generate_RockMass_Multi(dfn5)

generator.export_BlocksDDAOpt("DataSaveGeo_simp_juxtap1_3")
# for i in range(len(generator.blocks)):
# generator.blocks[i].export_BlockVtk("Block"+str(i))
