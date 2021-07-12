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
generator = Generator()

dfn = DFN()
dfn.set_FirstRecBlock([2, 2, 2], [4, 4, 3])
dfn.add_FractureSet()
generator.generate_RockMass_Multi(dfn)

# dfn = DFN()
# dfn.set_FirstRecBlock([6, 2, 1], [8, 4, 2])
# dfn.add_FractureSet()
# generator.generate_RockMass_Multi(dfn)

dfn = DFN()
dfn.set_FirstRecBlock([0, 0, 0], [10, 6, 1])
dfn.add_FractureSet()
generator.generate_RockMass_Multi(dfn)

generator.add_Fixed_Region([-.1, -.1, -.1], [20.1, 6.1, 1.1])
generator.export_BlocksDDAOpt("Test_freefall")
generator.export_BlocksVtk("Test_freefall")
