from numpy.core.shape_base import block
from numpy.lib.ufunclike import fix
import generateBlock as GB
import numpy as np
import copy

w = 4  # x
dx = w
l = 4  # y
dy = l
h = 2  # z
dz = h
# n = 101  # number of blocks, in one single dimension
# nd = 33
n = 14
nd = 4

blocks = []

vertexes = []
vertexes.append([0, l, h])
vertexes.append([0, 0, h])
vertexes.append([w, 0, h])
vertexes.append([w, l, h])
vertexes.append([0, l, 0])
vertexes.append([0, 0, 0])
vertexes.append([w, 0, 0])
vertexes.append([w, l, 0])
edges = [[0, 1], [1, 2], [2, 3], [3, 0], [0, 4], [1, 5], [2, 6], [3, 7],
         [4, 5], [5, 6], [6, 7], [7, 4]]
faces = [[0, 1, 2, 3], [0, 4, 5, 1], [0, 3, 7, 4], [2, 6, 7, 3],
         [1, 5, 6, 2], [4, 7, 6, 5]]

for iz in range(n):
    for ix in range(n):
        for iy in range(n):
            fixed = 0
            if iz == 0:
                fixed = 1
            vertexesi = copy.deepcopy(vertexes)
            if iz == 0 and ix > nd-1 and ix < n-nd and iy > nd-1 and iy < n-nd:
                continue
            if iz % 2 == 0:
                GB.move(vertexesi, [ix*dx, iy*dy, iz*dz])
                blocks.append(GB.Block(vertexesi, edges, faces, fixed=fixed))
            else:
                if ix < n-1 and iy < n-1:
                    GB.move(vertexesi, [ix*dx+0.5*dx, iy*dy+0.5*dy, iz*dz])
                    blocks.append(
                        GB.Block(vertexesi, edges, faces, fixed=fixed))

# GB.genBlk(blocks, "OneMillion.dda")
GB.genBlk(blocks, "multi_fall2k.dda")
