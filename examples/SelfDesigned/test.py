import generateBlock as GB
import numpy as np

b = 1
h = b*np.sqrt(3)/2
l = 5
vertexes = []
vertexes.append([0, 0.5*b, 0])
vertexes.append([0, 0, h])
vertexes.append([0, -0.5*b, 0])
vertexes.append([-l, 0, 0])
vertexes2 = [[0, l, -h], [-0.5*b, l, 0], [0.5*b, l, 0], [0, 0, 0]]
vertexes2 = GB.move(vertexes2, [-0.9*l, 0, (1-0.9)*h])

edges = [[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]]
faces = [[0, 1, 2], [0, 3, 1], [2, 1, 3], [2, 3, 0]]

tetra1 = GB.Block(vertexes, edges, faces, 1)
tetra2 = GB.Block(vertexes2, edges, faces, 0)


GB.genBlk([tetra1, tetra2], s="TwoSmallBlocks.dda")
