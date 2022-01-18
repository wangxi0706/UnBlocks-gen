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
off_x = 0
off_y = 0
off_z = 0
lx = 2
ly = 2
lz = 2
nz = 100

# dfn = DFN()
# dfn.set_RegionMaxCorner([lx*nx, ly, lz])
# dfn.add_FractureSet()
# for i in range(1, nx):
#     dfn.fractureSets[0].add_CircularFracture(
#         [lx*i, ly/2., lz/2.], 0, 90, lx)
generator = Generator()
for i in range(nz):
    dfn = DFN()
    min = [off_x+0, off_y+0, off_z+i*lz]
    max = [off_x+lx, off_y+ly, off_z+(i+1)*lz]
    dfn.set_FirstRecBlock(min, max)
    dfn.add_FractureSet()
    generator.input_RollerBound(i, [0,  1, 0])
    generator.input_RollerBound(i, [0,  -1, 0])
    # # for P wave propagation
    # generator.input_RollerBound(i, [0,  1, 0])
    # generator.input_RollerBound(i, [0,  -1, 0])
    # for S wave propagation
    generator.input_RollerBound(i, [-1, 0, 0])
    generator.input_RollerBound(i, [1, 0, 0])
    generator.generate_RockMass_Multi(dfn)


# for i in range(nz-1):
#     generator.input_BondedBlocks([i, i+1])

#################################################
# Boundary conditions
# be careful about the outer normal, should be outwards!
# generator.input_ViscousBound(0, [-1, 0, 0])
# generator.input_ViscousBound(nx-1, [1, 0, 0])
# generator.input_FixFaceBound(nz-1, [1, 0, 0])
# generator.input_InputBound(0, [-1, 0, 0])  # 0=disp, 1=load
# try input from the right
# generator.input_InputBound(nx-1, [1, 0, 0])  # 0=disp, 1=load
generator.add_Fixed_Region([-.1, -.1, -.1], [lx+0.1, ly+0.1, lz+0.1])
#################################################

# dfn.export_DFNVtk("dfnCreated")
# dfn.export_RegionVtk("modelRegion")

# generator = Generator()
# generator.generate_RockMass(dfn)
s = str(nz)+"PilingBlocks"
generator.export_BlocksVtk(s)
generator.export_BlocksDDAOpt(s)

# for i in range(len(generator.blocks)):
# generator.blocks[i].export_BlockVtk("Block"+str(i))
