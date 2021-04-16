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

# n=10, 1^2+...+n^2=(2n+1)n(n+1)/6

from unblocks import *

# dfn = DFN()
# n = 10
# dfn.set_RegionMaxMinCorner([8, 8, n], [0, 0, -n])
# dfn.add_FractureSet()
# generator = Generator()
# generator.add_Fixed_Region([-0.01, -0.01, -n-0.01], [10.1, 10.1, -n+1.01])
# dfn.fractureSets[0].add_CircularFracture([5, 5, 0], 0, 0, 1)
# for i in range(n-1):
#     dfn.fractureSets[0].add_CircularFracture([5, 5, i+1], 0, 0, 1)
#     dfn.fractureSets[0].add_CircularFracture([5, 5, -i-1], 0, 0, 1)
# generator.generate_RockMass(dfn)

from unblocks import *
n = 50
generator = Generator()
generator.add_Fixed_Region(
    [-5*n/2-.1, -5*n/2-.1, -.1], [5*n/2+.1, 5*n/2+.1, 2+.1])

dfns = []
for i in range(1, n+1):
    dfn = DFN()
    max = [5*i/2, 5*i/2, 2*(n-i+1)]
    min = [-5*i/2, -5*i/2, 2*(n-i)]
    dfn.set_RegionMaxMinCorner(max, min)
    print('max ', max, ' min ', min)
    dfn.add_FractureSet()
    for jj in range(0, i-1):
        j = jj - (i-2)/2
    # for j in range(int(-i/2+0.9), int(i/2+.1)):
        frac_cen_1 = [j*5., -5*i/2+.1, 2*(n-i+0.1)]
        frac_cen_2 = [j*5., 5*i/2-.1, 2*(n-i+0.1)]
        frac_cen_3 = [j*5., 5*i/2-.1, 2*(n-i+0.5)]
        print("frac point", frac_cen_1, frac_cen_2, frac_cen_3)
        dfn.fractureSets[0].add_TriangularFracture(
            frac_cen_1, frac_cen_2, frac_cen_3)
        frac_cen_1 = [-5*i/2+.1, j*5., 2*(n-i+0.1)]
        frac_cen_2 = [5*i/2-.1,  j*5., 2*(n-i+0.1)]
        frac_cen_3 = [5*i/2-.1,  j*5., 2*(n-i+0.5)]
        print("frac point", frac_cen_1, frac_cen_2, frac_cen_3)
        dfn.fractureSets[0].add_TriangularFracture(
            frac_cen_1, frac_cen_2, frac_cen_3)
    generator.generate_RockMass(dfn)
s = "Test_pyramid"+str(int((2*n+1)*n*(n+1)/6))

# s = "DataSaveGeo_BlkCut_20"
generator.export_BlocksDDAOpt(s)
generator.export_BlocksVtk(s)
