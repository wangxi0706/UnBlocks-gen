import generateBlock as GB
import numpy as np

h = 3

vertexes = []
vertexes.append([0, 0, 0])
vertexes.append([h, 0, 0])
vertexes.append([0, h, 0])
vertexes.append([0, 0, h])

edges = [[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]]
faces = [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]]

tetra1 = GB.Block(vertexes, edges, faces, 1)


GB.genBlk([tetra1], s="OneTetrahedron.dda")
