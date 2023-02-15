from unblocks import *
import generateBlock as GB
import numpy as np

################################
# For block B
w = 15  # x
l = 7.5  # y
h = 10  # z
t = 2.7

tol = 0
w += tol
l += tol
h += tol

vertexes = []
vertexes.append([tol, l-t+tol,h+tol])
vertexes.append([tol, tol,h+tol])
vertexes.append([w+tol, tol,h+tol])
vertexes.append([w+tol, l-t+tol,h+tol])
vertexes.append([tol, l+tol,tol])
vertexes.append([tol, tol,tol])
vertexes.append([w+tol, tol,tol])
vertexes.append([w+tol, l+tol,tol])
edges = [[0, 1], [1, 2], [2, 3], [3, 0], [0, 4], [1, 5], [2, 6], [3, 7],
         [4, 5], [5, 6], [6, 7], [7, 4]]
faces = [[0, 1, 2, 3], [0, 4, 5, 1], [0, 3, 7, 4], [2, 6, 7, 3],
         [1, 5, 6, 2], [4, 7, 6, 5]]
block1 = GB.Block(vertexes, edges, faces, fixed=1)

GB.genBlkTet([block1], "soil-small.poly")
