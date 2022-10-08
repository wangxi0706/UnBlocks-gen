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
nw = 1
nh = 10
nv = 10
w = 3
h = 3
v = 1.5
ww = nw*w
hh = nh*h
vv = nv*v
ww0 = ww/2-50
ww1 = ww/2+50
hh0 = hh/2-20
hh1 = hh/2+20
generator = Generator()

# two bases
nnh = 2.3
dfn = DFN()
min = [-1, -1, -v]
max = [ww+1, h*nnh, 0]
dfn.set_FirstRecBlock(min, max)
dfn.add_FractureSet()
generator.generate_RockMass_Multi(dfn)
dfn = DFN()
min = [-1, hh-h*nnh, -v]
max = [ww+1, hh+1, 0]
dfn.set_FirstRecBlock(min, max)
dfn.add_FractureSet()
generator.generate_RockMass_Multi(dfn)


for iw in range(nw):
    for iv in range(nv):
        if iv % 2 == 0:
            for ih in range(nh):
                dfn = DFN()
                min = [iw*w, ih*h, iv*v]
                max = [(iw+1)*w, (ih+1)*h, (iv+1)*v]
                dfn.set_FirstRecBlock(min, max)
                dfn.add_FractureSet()
                generator.generate_RockMass_Multi(dfn)
        else:
            dfn = DFN()
            min = [iw*w, (0)*h, iv*v]
            max = [(iw+1)*w, (.5)*h, (iv+1)*v]
            dfn.set_FirstRecBlock(min, max)
            dfn.add_FractureSet()
            generator.generate_RockMass_Multi(dfn)

            dfn = DFN()
            min = [iw*w, (nh-.5)*h, iv*v]
            max = [(iw+1)*w, (nh)*h, (iv+1)*v]
            dfn.set_FirstRecBlock(min, max)
            dfn.add_FractureSet()
            generator.generate_RockMass_Multi(dfn)
            for ih in range(nh-1):
                dfn = DFN()
                min = [iw*w, (ih+0.5)*h, iv*v]
                max = [(iw+1)*w, (ih+1.5)*h, (iv+1)*v]
                dfn.set_FirstRecBlock(min, max)
                dfn.add_FractureSet()
                generator.generate_RockMass_Multi(dfn)


generator.add_Fixed_Region([-1.1, -1.1, -v-0.1], [ww+1.1, hh+1.1, .1])
generator.export_BlocksDDAOpt("test_multi_fall"+str(nw*nh*nv))
generator.export_BlocksVtk("test_multi_fall"+str(nw*nh*nv))
