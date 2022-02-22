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

lx = 5
ly = 1
lz = 0.2

off_x = lx+0.5
off_y = ly+0.5

nXSlope=40
nYSlope=40

generator = Generator()
for ix in range(nXSlope):
    for iy in range(nYSlope):
        # bottom block
        dfn = DFN()
        min = [ix*off_x+ 0, iy*off_y+0, 0]
        max = [ix*off_x+ lx,iy*off_y+ ly, lz]
        dfn.set_FirstRecBlock(min, max)
        dfn.add_FractureSet()
        generator.generate_RockMass_Multi(dfn)
        # sliding block
        dfn = DFN()
        min = [ix*off_x+0.1*lx, iy*off_y+0.25*ly, lz]
        max = [ix*off_x+0.4*lx, iy*off_y+0.75*ly, lz+lz]
        dfn.set_FirstRecBlock(min, max)
        dfn.add_FractureSet()
        generator.generate_RockMass_Multi(dfn)

generator.add_Fixed_Region([-.1, -.1, -.1], [lx*2*nXSlope+0.1, ly*2*nYSlope+0.1, lz+0.1])
#################################################

s = str(nXSlope*nYSlope*2)+"MultiSlope"
generator.export_BlocksVtk(s)
generator.export_BlocksDDAOpt(s)

