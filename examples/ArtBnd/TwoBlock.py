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
lx = 2
ly = 2
lz = 2
nx = 2

# dfn = DFN()
# dfn.set_RegionMaxCorner([lx*nx, ly, lz])
# dfn.add_FractureSet()
# for i in range(1, nx):
#     dfn.fractureSets[0].add_CircularFracture(
#         [lx*i, ly/2., lz/2.], 0, 90, lx)
generator = Generator()
for i in range(nx):
    dfn = DFN()
    min = [i*lx, ly, lz]
    max = [(i+1)*lx, ly*2, lz*2]
    dfn.set_FirstRecBlock(min, max)
    dfn.add_FractureSet()
    generator.generate_RockMass_Multi(dfn)

# dfn.export_DFNVtk("dfnCreated")
# dfn.export_RegionVtk("modelRegion")

# generator = Generator()
# generator.generate_RockMass(dfn)
s = str(nx)+"blocks"
generator.export_BlocksVtk(s)
generator.export_BlocksDDAOpt(s)

# for i in range(len(generator.blocks)):
# generator.blocks[i].export_BlockVtk("Block"+str(i))
