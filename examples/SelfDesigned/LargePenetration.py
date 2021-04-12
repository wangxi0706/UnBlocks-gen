import generateBlock as GB
import numpy as np

################################
# For block B
w = 4  # x
l = 4  # y
h = 2  # z

vertexes = []
vertexes.append([-w/2, l/2, 0])
vertexes.append([-w/2, -l/2, 0])
vertexes.append([w/2, -l/2, 0])
vertexes.append([w/2, l/2, 0])
vertexes.append([-w/2, l/2, -h])
vertexes.append([-w/2, -l/2, -h])
vertexes.append([w/2, -l/2, -h])
vertexes.append([w/2, l/2, -h])
edges = [[0, 1], [1, 2], [2, 3], [3, 0], [0, 4], [1, 5], [2, 6], [3, 7],
         [4, 5], [5, 6], [6, 7], [7, 4]]
faces = [[0, 1, 2, 3], [0, 4, 5, 1], [0, 3, 7, 4], [2, 6, 7, 3],
         [1, 5, 6, 2], [4, 7, 6, 5]]
block1 = GB.Block(vertexes, edges, faces, fixed=1)


#######################################
# For block A
w = w/2
l = l/2
h = h/2

vertexes = []
vertexes.append([-w/2, l/2, 0])
vertexes.append([-w/2, -l/2, 0])
vertexes.append([w/2, -l/2, 0])
vertexes.append([w/2, l/2, 0])
vertexes.append([-w/2, l/2, -h])
vertexes.append([-w/2, -l/2, -h])
vertexes.append([w/2, -l/2, -h])
vertexes.append([w/2, l/2, -h])
GB.move(vertexes, [0, 0, h-0.05])

block2 = GB.Block(vertexes, edges, faces)

#######################################
GB.genBlk([block1, block2], "LargePenetration.dda")
